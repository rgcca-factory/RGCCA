#' @importFrom MASS ginv

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
  x$h_tilde <- x$QMX %*% grad
  mu_max <- 0.5 * sqrt(sum(x$h_tilde**2))
  mu_min <- max(0, sapply(seq_len(x$p), function(k) {0.5 * sqrt(sum(x$h_tilde[k:x$p]**2)) - max(x$D[k:x$p])}))
  L <- function(mu) {0.5 * sum(x$h_tilde**2 / (x$D + 2 * mu)) + mu}
  x$mu <- optimize(f = L, interval = c(mu_min, mu_max), tol = .Machine$double.eps)$minimum
  
  x$a <- - x$MQ %*% (x$h_tilde / (x$D + 2 * x$mu))
  return(block_project(x))
}

#' @export
block_update.dual_ac_block <- function(x, grad) {
  x$h <- x$h_K %*% grad
  mu_max <- 1E10 #TODO find a way to define an upper bound for mu in the dual setting (and a lower bound?)
  L <- function(mu) {
    0.5 * drop(t(x$h) %*% ginv(x$B + 2 * mu * x$KM, tol = .Machine$double.eps * x$n) %*% x$h) + mu
  }
  x$mu <- optimize(f = L, interval = c(0, mu_max), tol = .Machine$double.eps)$minimum

  x$alpha <- - ginv(x$B + 2 * x$mu * x$KM, tol = .Machine$double.eps * x$n) %*% x$h
  return(block_project(x))
}