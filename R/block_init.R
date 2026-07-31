#' @importFrom MASS ginv
#' @importFrom RSpectra eigs_sym

block_init <- function(x, init = "svd") {
  UseMethod("block_init")
}

#' @export
block_init.block <- function(x, init = "svd") {
  if (init == "svd") {
    x$a <- initsvd(x$x, dual = FALSE)
  } else {
    x$a <- rnorm(x$p)
  }
  
  return(block_project(x))
}

#' @export
block_init.dual_block <- function(x, init = "svd") {
  if (init == "svd") {
    x$alpha <- initsvd(x$x, dual = TRUE)
  } else {
    x$alpha <- rnorm(x$n)
  }
  
  return(block_project(x))
}

#' @export
block_init.primal_regularized_block <- function(x, init = "svd") {
  x$M <- ginv(
    x$tau * diag(x$p) + (1 - x$tau) * pm(t(x$x), x$x, na.rm = x$na.rm) / x$N
  )
  NextMethod()
}

#' @export
block_init.dual_regularized_block <- function(x, init = "svd") {
  x$M <- ginv(x$tau * diag(x$n) + (1 - x$tau) * x$K / x$N)
  NextMethod()
}

#' @export
block_init.ac_block <- function(x, init = "svd") {
  x$M <- x$tau * diag(x$p) + (1 - x$tau) * pm(t(x$x), x$x, na.rm = x$na.rm) / x$N
  if (x$tau == 1) {
    x$sqrt_M <- diag(x$p)
    x$sqrt_M_inv <- diag(x$p)
  } else {
    x$sqrt_M <- sqrt_matrix(x$M)
    x$sqrt_M_inv <- sqrt_matrix(x$M, inv = T)
  }
  
  P <- pm(x$x, x$sqrt_M_inv, na.rm = x$na.rm)
  x$B <- 2 * 1/x$N * x$gamma_confounders * pm(
    t(P), pm(
      x$confounders, P, na.rm = x$na.rm), 
    na.rm = x$na.rm)
  
  eigen_dec <- eigen(x$B, symmetric = T)
  x$D <- eigen_dec$values
  x$Q <- eigen_dec$vectors
  
  x$QMX <- - 2/x$N * t(x$Q) %*% t(P) 
  
  x$MQ <- x$sqrt_M_inv %*% x$Q
  
  NextMethod()
}

#' @export
block_init.dual_ac_block <- function(x, init = "svd") {
  x$M_n <- x$tau * diag(x$n) + (1 - x$tau)/x$N * x$K 
  x$KM <- pm(x$K, x$M_n, na.rm = x$na.rm)
  
  x$B <- 2/x$N * x$gamma_confounders * x$K %*% x$confounders %*% x$K
  x$h_K <- -2/x$N * x$K
  
  NextMethod()  
}
