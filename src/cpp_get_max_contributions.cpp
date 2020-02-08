#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <cstdbool>
#include <cstddef>
#include <numeric>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// compare-function
bool compare(std::vector<double> const &rowA, std::vector<double> const &rowB) {
  // If the values of the second-last column are not equal,
  // use them to order rowA and rowB.
  int varIndex1 = rowA.size() - 2;
  if (rowA[varIndex1] != rowB[varIndex1]) {
    return rowA[varIndex1] > rowB[varIndex1];
  } else {
    // If the values of the second-last column are equal,
    // use the values of the first column to order rowA and rowB.
    return rowA[0] > rowB[0];
  }
}

// results container
struct resList {
  std::string varName;

  std::vector<int> uw_ids;
  std::vector<double> uw_vals;
  double uw_spread;
  double uw_sum;
  double uw_mean;

  std::vector<int> w_ids;
  std::vector<double> w_vals;
  double w_spread;
  double w_sum;
  double w_mean;

  bool even_contributors;
};

// [[Rcpp::export]]
List cpp_get_max_contributions(List logicals_R, NumericMatrix microdat_R,
                               String wvar_R, CharacterVector nv_R,
                               int top_k_in, int n_threads) {

  std::size_t nLogicals = logicals_R.size();
  std::size_t nNumVar = nv_R.size();
  std::size_t nRow = microdat_R.nrow();
  std::size_t nCol = microdat_R.ncol();

  // create C++ objects from R inputs
  std::vector<std::vector<double>> logicals(nLogicals);
  for (std::size_t i = 0; i < nLogicals; i++) {
    logicals[i] = as<std::vector<double>>(logicals_R[i]);
  }
  std::vector<std::vector<double>> microdat(nRow, std::vector<double>(nCol));
  for (std::size_t i = 0; i < nRow; i++) {
    for (std::size_t j = 0; j < nCol; j++) {
      microdat[i][j] = microdat_R(i, j);
    }
  }
  std::string wvar = wvar_R;
  std::vector<std::string> nv = as<std::vector<std::string>>(nv_R);

  // create list for each set of indices
  std::vector<std::vector<resList>> res(nLogicals,
                                        std::vector<resList>(nNumVar));

// for each set of indices:
#ifdef _OPENMP
#pragma omp parallel for num_threads(n_threads)
#endif
  for (std::size_t i = 0; i < nLogicals; i++) {
    // create subset of data, subsetted by indices
    std::vector<double> keep = logicals[i];
    std::size_t nKeep = std::accumulate(keep.begin(), keep.end(), 0);
    // ordering of variables:
    // 0: weight-variable (defined in wvar),
    // 1: .tmpid,
    // 2 to (nNumVar + 1): numeric variables specified in "nv"
    std::vector<std::vector<double>> XX(nKeep, std::vector<double>(nCol + 2));

    std::size_t l = 0;
    for (std::size_t j = 0; j < (std::size_t)keep.size(); j++) {
      if (keep[j] == 1) {
        for (std::size_t k = 0; k < nCol; k++) {
          XX[l][k] = microdat[j][k];
        }
        l = l + 1;
      }
    }

    int top_k = std::min(top_k_in, (int)XX.size());

    // create variable, sort data, compute contributions
    for (std::size_t v = 0; v < nNumVar; v++) {
      // create new variables
      for (std::size_t k = 0; k < nKeep; k++) {
        XX[k][1 + nNumVar + 1] = std::abs(XX[k][v + 2]);
        XX[k][1 + nNumVar + 2] = XX[k][v + 2] * XX[k][0];
      }
      // original R implementation reorders matrix in-place (doing the same here
      // for comparability). reordering based on "fresh" matrices for each
      // numvar yields different results.
      // sort data by newly created variables
      std::stable_sort(XX.begin(), XX.end(), compare);

      // compute contributions etc.
      std::vector<int> uw_ids(top_k);
      std::vector<double> uw_vals(top_k);
      double uw_spread;
      double uw_sum;
      double uw_mean;
      bool even_contributors;

      std::vector<int> w_ids = uw_ids;
      std::vector<double> w_vals(top_k);
      double w_spread;
      double w_sum;
      double w_mean;

      if (top_k == 0) {
        uw_ids.resize(1);
        uw_vals.resize(1);
        w_ids.resize(1);
        w_vals.resize(1);

        uw_ids[0] = XX[0][1];
        uw_vals[0] = 0;
        uw_spread = 0;
        uw_sum = 0;
        uw_mean = 0;

        w_ids[0] = uw_ids[0];
        w_vals[0] = uw_vals[0];
        w_spread = 0;
        w_sum = 0;
        w_mean = 0;
        even_contributors = 0;
      } else {
        // unweighted
        for (int j = 0; j < top_k; j++) {
          uw_ids[j] = XX[j][1];
          uw_vals[j] = XX[j][v + 2];
        }
        std::vector<double> helper(XX.size());
        for (std::size_t j = 0; j < XX.size(); j++) {
          helper[j] = XX[j][v + 2];
        }
        double uw_min = *std::min_element(helper.begin(), helper.end());
        double uw_max = *std::max_element(helper.begin(), helper.end());
        uw_spread = uw_max - uw_min;
        uw_sum = std::accumulate(helper.begin(), helper.end(), 0);
        uw_mean = uw_sum / XX.size();

        // weighted
        for (int j = 0; j < top_k; j++) {
          w_vals[j] = XX[j][1 + nNumVar + 2];
          w_ids[j] = uw_ids[j];
        }
        std::vector<double> wvar(XX.size());
        for (std::size_t j = 0; j < XX.size(); j++) {
          helper[j] = XX[j][1 + nNumVar + 2];
          wvar[j] = XX[j][0];
        }
        double w_min = *std::min_element(helper.begin(), helper.end());
        double w_max = *std::max_element(helper.begin(), helper.end());
        w_spread = w_max - w_min;
        w_sum = std::accumulate(helper.begin(), helper.end(), 0);
        w_mean = w_sum / (double)std::accumulate(wvar.begin(), wvar.end(), 0);

        // we compute if the number of contributors to the cell
        // is even or odd. This information can later be used if
        // we have different ptables (parity-case)
        even_contributors = (XX.size() % 2 == 0);
      }
      // results list
      res[i][v].varName = nv[v];
      res[i][v].uw_ids = uw_ids;
      res[i][v].w_ids = w_ids;
      res[i][v].uw_spread = uw_spread;
      res[i][v].uw_sum = uw_sum;
      res[i][v].uw_mean = uw_mean;
      res[i][v].w_spread = w_spread;
      res[i][v].w_sum = w_sum;
      res[i][v].w_mean = w_mean;
      res[i][v].uw_vals = uw_vals;
      res[i][v].w_vals = w_vals;
      res[i][v].even_contributors = even_contributors;
    }
  }

  // create results object for R
  List res_R(nLogicals);
  res_R.names() = logicals_R.names();
  for (std::size_t i = 0; i < nLogicals; i++) {
    List out(nNumVar);
    out.names() = nv_R;
    for (std::size_t v = 0; v < nNumVar; v++) {
      out[v] = List::create(
          _["uw_ids"] = res[i][v].uw_ids, _["w_ids"] = res[i][v].w_ids,
          _["uw_spread"] = res[i][v].uw_spread, _["uw_sum"] = res[i][v].uw_sum,
          _["uw_mean"] = res[i][v].uw_mean, _["w_spread"] = res[i][v].w_spread,
          _["w_sum"] = res[i][v].w_sum, _["w_mean"] = res[i][v].w_mean,
          _["uw_vals"] = res[i][v].uw_vals, _["w_vals"] = res[i][v].w_vals,
          _["even_contributors"] = res[i][v].even_contributors);
    }
    res_R[i] = out;
  }
  return res_R;
}
