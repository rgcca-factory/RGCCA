format_bootstrap_list <- function(W, rgcca_res, n_boot, n = 1) {
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
  pjs <- vapply(rgcca_res$call$blocks, ncol, FUN.VALUE = integer(1))
  list_res_X <- array(unlist(lapply(X, function(x) Reduce(rbind, x))),
    dim = c(sum(pjs), ndefl_max, n_boot)
  )
  rownames(list_res_X) <- unlist(lapply(
    rgcca_res$call$blocks, colnames
  ))
  f <- unlist(mapply(function(x, y) rep(x, each = y), seq(pjs), pjs))
  list_res_X <- lapply(seq(ndefl_max), function(i) {
    x <- split(list_res_X[, i, ], f)
    x <- lapply(seq_along(x), function(j) {
      y <- matrix(x[[j]], nrow = pjs[j])
      rownames(y) <- colnames(rgcca_res$call$blocks[[j]])
      return(y)
    })
    names(x) <- names(rgcca_res$call$blocks)
    return(x)
  })
}
