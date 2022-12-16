#' Parallel lapply with progress bar
#' @inheritParams rgcca_cv
#' @param X a vector (atomic or list) or an expression object.
#' @param FUN the function to be applied to each element of X.
#' @param ... optional arguments to FUN.
#' @return The result of lapply(X, FUN, ...)
#' @noRd
par_pblapply <- function(X, FUN, ..., n_cores = 1, verbose = TRUE) {
  check_integer("n_cores", n_cores, min = 0)

  if (!verbose) {
    pbapply::pboptions(type = "none")
  } else {
    pbapply::pboptions(type = "timer")
  }

  is_windows <- Sys.info()["sysname"] == "Windows"

  if (n_cores <= 1) {
    cl <- NULL
  } else if (is_windows) {
    cl <- parallel::makeCluster(n_cores)
    parallel::clusterExport(cl, NULL, envir = environment())
  } else {
    cl <- n_cores
  }

  W <- pbapply::pblapply(X, FUN, ..., cl = cl)

  if (is_windows && !is.null(cl)) parallel::stopCluster(cl)

  return(W)
}
