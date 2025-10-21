block_postprocess <- function(x, ctrl) {
  UseMethod("block_postprocess")
}

#' @export
block_postprocess.block <- function(x, ctrl) {
  if (ctrl && (x$a[1] < 0)) {
    x$a <- -x$a
    x$Y <- -x$Y
  }
  return(x)
}

#' @export
block_postprocess.sparse_block <- function(x, ctrl) {
  l2_sat <- norm(x$a, "2")
  if (abs(l2_sat - 1) > x$tol) {
    if (l2_sat < .Machine$double.eps) {
      warning(
        "The l2 norm of the block weight vector #",
        x$j, " is too small :", l2_sat
      )
    } else {
      warning(
        "The l2 constraint is not saturated for block #", x$j,
        ". The intersection of the l1 and l2 spheres is empty for ",
        "a sparsity parameter equal to ", x$sparsity,
        ". Try to increase the value of the sparsity parameter."
      )
    }
  }
  NextMethod()
}

#' @export
block_postprocess.separable_regularized_tensor_block <- function(x, ctrl) {
  x$factors <- lapply(seq_along(x$factors), function(m) {
    x$M[[m]] %*% x$factors[[m]]
  })
  x$a <- Reduce(khatri_rao, rev(x$factors)) %*% x$lambda
  NextMethod()
}