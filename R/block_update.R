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
    #mu_max <- 1/2 * sum(x$e**2 / x$d)
    #x$mu <- mean(gtools::binsearch(fun = function(mu) {sum(x$e**2 / (x$d + 2 * mu)**2) - 1},range = c(0, mu_max))$where)
    
    L <- function(mu) {1/2 * sum(x$e**2 / (x$d + 2 * mu)) + mu}
    grad_L <- function(mu) {- sum(x$e**2 / (x$d + 2 * mu)**2) + 1}
    res <- optim(par = 0, L, grad_L, method = "BFGS")
    x$mu <- res$par
  } else if (x$algo == 3) {
    x$h_tilde <- pm(x$h_tilde_QMX, grad, na.rm = x$na.rm)
    mu_max <- 0.5 * sum(x$h_tilde**2)
    L <- function(mu) {0.5 * sum(x$h_tilde**2 / (x$eigen_val_Bplus1 + 2 * mu - 1)) + mu}
    x$mu <- optimize(f = L, interval = c(0, mu_max))$minimum
  }
  return(block_project(x))
}