#' Get the rgcca arguments from a fitted object and default arguments.
#' @noRd
get_rgcca_args <- function(object, default_args = list()) {
  if (class(object) %in% c("rgcca", "permutation", "cval")) {
    rgcca_args <- object$call

    if (class(object) %in% c("permutation", "cval")) {
      if (object$par_type == "tau") {
        rgcca_args$tau <- object$bestpenalties
      }
      if (object$par_type == "ncomp") {
        rgcca_args$ncomp <- object$bestpenalties
      }
      if (object$par_type == "sparsity") {
        rgcca_args$sparsity <- object$bestpenalties
      }
    }
  } else {
    rgcca_args <- list(
      tau = default_args$tau,
      tol = default_args$tol,
      init = default_args$init,
      bias = default_args$bias,
      quiet = default_args$quiet,
      scale = default_args$scale,
      ncomp = default_args$ncomp,
      blocks = default_args$blocks,
      scheme = default_args$scheme,
      method = default_args$method,
      verbose = default_args$verbose,
      sparsity = default_args$sparsity,
      response = default_args$response,
      NA_method = default_args$NA_method,
      n_iter_max = default_args$n_iter_max,
      connection = default_args$connection,
      superblock = default_args$superblock,
      scale_block = default_args$scale_block
    )
  }
  return(rgcca_args)
}
