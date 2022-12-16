#' Format bootstrap results
#'
#' Raw bootstrap results are stored in a list of n_boot elements where each
#' element of this list is a list of J matrices. Each of these matrices has
#' ncomp[j] columns where ncomp[j] is the number of components for block j.
#' The aim of this function is to reorganize this list to have a list of
#' max(ncomp) elements, where each each element is a list of J matrices of
#' n_boot columns.
#' @param W raw bootstrap results
#' @param rgcca_res a fitted rgcca object
#' @param n_boot the number of bootstrap samples
#' @param n integer (1 or 2) indicating which of raw bootstrap list to reorder
#' @return Reordered list of bootstrap results.
#' @noRd
format_bootstrap_list <- function(W, rgcca_res, n_boot, n = 1) {
  J <- length(rgcca_res$call$blocks)
  ndefl_max <- max(rgcca_res$call$ncomp)
  X <- lapply(W, `[[`, n)
  # Add columns of NA for missing components
  X <- lapply(X, function(x) {
    lapply(seq_along(x), function(j) {
      cbind(
        x[[j]],
        matrix(NA, nrow(x[[j]]), ndefl_max - rgcca_res$call$ncomp[j])
      )
    })
  })
  # Change order of elements
  pjs <- vapply(rgcca_res$blocks, ncol, FUN.VALUE = integer(1))
  list_res_X <- array(unlist(lapply(X, function(x) Reduce(rbind, x))),
    dim = c(sum(pjs), ndefl_max, n_boot)
  )
  rownames(list_res_X) <- unlist(lapply(
    rgcca_res$blocks, colnames
  ))
  f <- unlist(Map(function(x, y) rep(x, each = y), seq(pjs), pjs))
  list_res_X <- lapply(seq(ndefl_max), function(i) {
    # Split back into blocks and remove superblock
    x <- split(list_res_X[, i, ], f)[-(J + 1)]
    x <- lapply(seq_along(x), function(j) {
      y <- matrix(x[[j]], nrow = pjs[j])
      rownames(y) <- colnames(rgcca_res$blocks[[j]])
      return(y)
    })
    names(x) <- names(rgcca_res$call$blocks)
    return(x)
  })
}
