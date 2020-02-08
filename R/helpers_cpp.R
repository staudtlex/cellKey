# for each numerical variable and each cell, get the top_k contributions
# along with its ids and its values
#' @importFrom parallel detectCores
#' @importFrom Rcpp sourceCpp
#' @useDynLib cellKey, .registration = TRUE
.cpp_get_max_contributions <- function(indices, microdat, wvar, nv, top_k, cores = 1) {
  varlist <- c(wvar, ".tmpid", nv)
  microdat <- as.matrix(microdat[, varlist, with = FALSE])
  keep <- lapply(indices, function(x) microdat[, ".tmpid"] %in% x)
  if (cores > detectCores()) {
    cores = detectCores()
    warning("Too many cores requested. Resetting to number of available cores.\n")
  } else if (cores < 1) {
    stop("Number of cores must be larger than 0.\n")
  }
  cpp_get_max_contributions(
    logicals_R = keep, 
    microdat_R = microdat, 
    wvar_R = wvar, 
    nv_R = nv, 
    top_k_in = top_k,
    n_threads = cores)
}

.get_max_contributions <- function(indices, microdat, wvar, nv, top_k) {
  res <- vector("list", length = length(indices))
  names(res) <- names(indices)
  
  for (i in seq_len(length(res))) {
    out <- vector("list", length = length(nv))
    names(out) <- nv
    xx <- microdat[microdat$.tmpid %in% indices[[i]]]
    top_k <- min(top_k, nrow(xx))
    for (v in nv) {
      xx$.tmpordervar <- abs(xx[[v]])
      xx$.tmpweightvar <- xx[[v]] * xx[[wvar]]
      setorderv(xx, c(".tmpordervar", wvar), order = c(-1L, -1L))
      
      # unweighted
      out[[v]]$uw_ids <- xx$.tmpid[1:top_k]
      out[[v]]$w_ids <- out[[v]]$uw_ids
      
      if (top_k == 0) {
        out[[v]]$uw_spread <- out[[v]]$w_spread <- 0
        out[[v]]$uw_sum <- out[[v]]$w_sum <- 0
        out[[v]]$uw_mean <- out[[v]]$w_mean <- 0
        out[[v]]$uw_vals <- out[[v]]$w_vals <- 0
      } else {
        out[[v]]$uw_spread <- diff(range(xx[[v]], na.rm = TRUE))
        out[[v]]$uw_sum <- sum(xx[[v]], na.rm = TRUE)
        out[[v]]$uw_mean <- out[[v]]$uw_sum / nrow(xx)
        
        out[[v]]$w_spread <- diff(range(xx$.tmpweightvar, na.rm = TRUE))
        out[[v]]$w_sum <- sum(xx$.tmpweightvar, na.rm = TRUE)
        out[[v]]$w_mean <- out[[v]]$w_sum / sum(xx[[wvar]], na.rm = TRUE)
        
        out[[v]]$uw_vals <- xx[[v]][1:top_k]
        out[[v]]$w_vals <- xx$.tmpweightvar[1:top_k]
      }
      # we compute if the number of contributors to the cell
      # is even or odd. This information can later be used if
      # we have different ptables (parity-case)
      out[[v]]$even_contributors <- nrow(xx) %% 2 == 0
    }
    res[[i]] <- out
  }
  res
}