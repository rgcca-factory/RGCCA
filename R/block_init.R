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
  if (x$algo == 1) {
    x$M <- x$tau * diag(x$p) + (1 - x$tau) * pm(t(x$x), x$x, na.rm = x$na.rm) / x$N
    x$M_inv <- ginv(x$M)
    x$sqrt_M <- sqrt_matrix(x$M)
    x$sqrt_M_inv <- sqrt_matrix(x$M_inv) #TODO check if result is the same with sqrt_matrix(x$M, inv = T)
  
    P <- pm(x$x, x$sqrt_M_inv, na.rm = x$na.rm)
    x$B <- pm(
      t(P), pm(
        x$confounders, P, na.rm = x$na.rm), 
      na.rm = x$na.rm)
    
    x$mu <- x$penalty_coef *
      RSpectra::eigs_sym(A = x$B, k = 1, which = "LA", opts = list(retvec = F))$values
    
    x$f_left <- 1/(2 * x$mu) * t(P)
    
    x$f_right <- x$penalty_coef / x$mu * pm(x$B, x$sqrt_M, na.rm = x$na.rm)
  } else if (x$algo == 2) {
    x$M <- x$tau * diag(x$p) + (1 - x$tau) * pm(t(x$x), x$x, na.rm = x$na.rm) / x$N
    x$M_inv <- ginv(x$M)
    x$sqrt_M <- sqrt_matrix(x$M)
    x$sqrt_M_inv <- sqrt_matrix(x$M_inv)
    
    P <- pm(x$x, x$sqrt_M_inv, na.rm = x$na.rm)
    
    res_svd <- svd(2 * x$penalty_coef * pm(
      t(P), pm(
        x$confounders, P, na.rm = x$na.rm), 
      na.rm = x$na.rm))
    x$d <- res_svd$d
    x$Q <- res_svd$v
    
    x$e_QM <- - pm(x$Q, x$sqrt_M, na.rm = x$na.rm)
    x$a_MQ <- - pm(x$sqrt_M_inv, t(x$Q), na.rm = x$na.rm)
  }
  
  NextMethod()
}
