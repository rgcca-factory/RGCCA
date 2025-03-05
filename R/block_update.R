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
block_update.ac_block <- function(x, grad) {
  x$f <- pm(x$f_left, 1/x$N * grad, na.rm = x$na.rm) -
    pm(x$f_right, x$a, na.rm = x$na.rm)
  return(block_project(x))
}