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
block_update.graphnet_block <- function(x, grad) {
  x$a <- as.matrix(
    pm(t(x$x), grad, na.rm = x$na.rm)  +
    x$lambda * x$graph_laplacian %*% x$a
  ) # if graph_laplacian is a sparse matrix
  # it a necessary to cast the result as a dense matrix
  return(block_project(x))
}
