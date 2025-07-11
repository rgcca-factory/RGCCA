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
  if (x$algo == 3) {
    x$M <- x$tau * diag(x$p) + (1 - x$tau) * pm(t(x$x), x$x, na.rm = x$na.rm) / x$N
    if (x$tau == 1) {
      x$sqrt_M <- diag(x$p)
      x$sqrt_M_inv <- diag(x$p)
    } else {
      x$sqrt_M <- sqrt_matrix(x$M)
      x$sqrt_M_inv <- sqrt_matrix(x$M, inv = T)
    }
    
    P <- pm(x$x, x$sqrt_M_inv, na.rm = x$na.rm)
    x$B <- 2 * 1/x$N * x$penalty_coef * pm(
      t(P), pm(
        x$confounders, P, na.rm = x$na.rm), 
      na.rm = x$na.rm)
    
    eigen_dec <- eigen(x$B, symmetric = T)
    x$D <- eigen_dec$values
    x$Q <- eigen_dec$vectors
    
    # cat("smallest eigenval= ", min(eigen_dec$values), "\n")
    # 
    # x$rank_B <- sum(x$D > max(x$D) * .Machine$double.eps * x$p)
    # 
    x$QMX <- - 2/x$N * t(x$Q) %*% t(P) 
    
    x$MQ <- x$sqrt_M_inv %*% x$Q
  } else if (x$algo == 5) {
    if (x$tau == 1) {
      x$M <- diag(x$p)
      x$sqrt_M <- diag(x$p)
      x$sqrt_M_inv <- diag(x$p)
    } else {
      x$M <- x$tau * diag(x$p) + (1 - x$tau) * pm(t(x$x), x$x, na.rm = x$na.rm) / x$N
      x$sqrt_M <- sqrt_matrix(x$M)
      x$sqrt_M_inv <- sqrt_matrix(x$M, inv = T)
    }
    
    P <- pm(x$x, x$sqrt_M_inv, na.rm = x$na.rm)
    x$B <- 2 * 1/x$N * x$penalty_coef * pm(
      t(P), pm(
        x$confounders, P, na.rm = x$na.rm), 
      na.rm = x$na.rm)
    
    eigen_dec <- eigen(x$B, symmetric = T)
    rank <- sum(eigen_dec$values > x$p * max(eigen_dec$values) * .Machine$double.eps)
    x$D <- c(eigen_dec$values[1:rank], rep(0, x$p - rank))
    x$Q <- eigen_dec$vectors
    
    x$MQ <- x$sqrt_M_inv %*% x$Q
    x$QMX <- - 2/x$N * t(x$Q) %*% t(P) 
    x$gamma <- 0.9 * 2 / x$D[1]
  }
  
  NextMethod()
}

#' @export
block_init.dual_ac_block <- function(x, init = "svd") {
  if (x$algo == 3) {
    x$M_n <- x$tau * diag(x$n) + (1 - x$tau)/x$N * x$K 
    #x$M_n_inv <- ginv(x$M_n)
    x$KM <- pm(x$K, x$M_n, na.rm = x$na.rm)
    
    x$B <- 2/x$N * x$penalty_coef * x$K %*% x$confounders %*% x$K
    x$h_K <- -2/x$N * x$K

  } else if (x$algo == 5) {
    x$M_n <- x$tau * diag(x$n) + (1 - x$tau)/x$N * x$K 
    x$KM <- pm(x$K, x$M_n, na.rm = x$na.rm)
    
    x$B <- 2/x$N * x$penalty_coef * x$K %*% x$confounders %*% x$K
    eigen_dec_B <- eigen(x$B, symmetric = T)
    rank_B <- sum(eigen_dec_B$values > x$n * max(eigen_dec_B$values) * .Machine$double.eps)
    x$D_B <- c(eigen_dec_B$values[1:rank_B], rep(0, x$n - rank_B))
    x$Q_B <- eigen_dec_B$vectors
    
    x$A <- t(x$Q) %*% x$KM %*% x$Q
    eigen_dec_A <- eigen(x$A, symmetric = T)
    rank_A <- sum(eigen_dec_A$values > x$n * max(eigen_dec_A$values) * .Machine$double.eps)
    x$D_A <- c(eigen_dec_A$values[1:rank_A], rep(0, x$n - rank_A))
    x$Q_A <- eigen_dec_A$vectors
    
    x$QK <- -2/x$N * t(x$Q_B) %*% x$K
    x$gamma <- 0.9 * 2 / x$D_B[1]
    x$QKMQ <- t(x$Q_B) %*% x$KM %*% x$Q_B
    
  }
  NextMethod()  
}
