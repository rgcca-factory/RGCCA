block_update <- function(x, grad) {
  UseMethod("block_update")
}

#' @export
block_update.block <- function(x, grad) {
  x$a <- pm(t(x$x), grad, na.rm = x$na.rm)
  return(block_project(x))
}

#' @export
block_update.dual_block <- function(x, grad) {
  x$alpha <- grad
  return(block_project(x))
}

#' @export
block_update.tensor_block <- function(x, grad) {
  grad <- array(
    pm(t(matrix(x$x, nrow = nrow(x$x))), grad, na.rm = x$na.rm),
    dim = dim(x$x)[-1]
  )
  other_factors <- NULL
  # Update factors
  for (m in seq_along(dim(x$x)[-1])) {
    grad_m <- matrix(
      aperm(grad, c(m, seq_along(dim(grad))[-m])), nrow = dim(x$x)[m + 1]
    )
    grad_m <- grad_m %*% khatri_rao(
      Reduce(khatri_rao, rev(x$factors[-seq_len(m)])), other_factors
    )
    x$factors[[m]] <- grad_m %*% diag(x$weights, nrow = length(x$weights))
    x$factors[[m]] <- x$factors[[m]] / norm(x$factors[[m]], type = "2")

    other_factors <- khatri_rao(x$factors[[m]], other_factors)
  }
  # Update weights
  x$weights <- t(other_factors) %*% as.vector(grad)
  x$weights <- drop(x$weights) / norm(drop(x$weights), type = "2")
  return(block_project(x))
}
