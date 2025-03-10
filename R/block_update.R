#' @importFrom gtools binsearch

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
  if (x$algo == 1) {
    x$f <- pm(x$f_left, 1/x$N * grad, na.rm = x$na.rm) -
      pm(x$f_right, x$a, na.rm = x$na.rm)
  } else if (x$algo == 2) {
    x$e <- x$e_QM %*% pm(t(x$x), 1/x$N * grad, na.rm = x$na.rm)
    mu_max <- 1/2 * drop(sweep(t(x$e), 2, (x$d + 1E-25)**(-1), "*") %*% x$e)
    x$mu <- mean(gtools::binsearch(
      fun = function(mu) {drop(sweep(t(x$e), 2, (x$d + 2 * mu)**(-2), "*") %*% x$e) - 1},
      range = c(0, mu_max))$where)
  }
  return(block_project(x))
}