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
    #x$M_inv <- ginv(x$M)
    x$sqrt_M <- sqrt_matrix(x$M)
    x$sqrt_M_inv <- sqrt_matrix(x$M, inv = T)
  
    P <- pm(x$x, x$sqrt_M_inv, na.rm = x$na.rm)
    x$B <- 1/x$N * pm(
      t(P), pm(
        x$confounders, P, na.rm = x$na.rm), 
      na.rm = x$na.rm)
    
    x$mu <- x$penalty_coef *
      RSpectra::eigs_sym(A = x$B, k = 1, which = "LA", opts = list(retvec = F))$values
    
    x$f_left <- 1/(2 * x$mu) * t(P)
    
    x$f_right <- x$penalty_coef / x$mu * pm(x$B, x$sqrt_M, na.rm = x$na.rm)
  } else if (x$algo == 3) {
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
    
    cat("smallest eigenval= ", min(eigen_dec$values), "\n")
    
    x$rank_B <- sum(x$D > max(x$D) * .Machine$double.eps * x$p)
    
    x$h_tilde_QMX <- - 1/x$N * pm(
      t(x$Q),
      pm(x$sqrt_M_inv, 
         t(x$x), na.rm = x$na.rm), 
      na.rm = x$na.rm)
    
    x$a_MQ <- - pm(x$sqrt_M_inv, x$Q, na.rm = x$na.rm)
  } else if (x$algo == 4) {
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
    
    cat("smallest eigenval= ", min(eigen_dec$values), "\n")
    
    # Check that the matrix is positive
    #if (any(eigen_dec$values < 0)) {
    #  D <- eigen_dec$values
    #  D_tilde <- D[D > x$p * max(D) * .Machine$double.eps]
    #  #Q_tilde <- eigen_dec$vectors[, 1:length(D)]
    #  
    #  D_tilde_complete <- c(D_tilde, rep(1E-7, x$p - length(D_tilde)))
    #  x$B <- eigen_dec$vectors %*% diag(D_tilde_complete, nrow = x$p) %*% t(eigen_dec$vectors)
    #}
    
    #res_svd <- svd(x$B)
    #svd_D_tilde <- res_svd$d[res_svd$d > x$p * max(res_svd$d) * .Machine$double.eps]
    #u_tilde <- res_svd$u[, 1:length(svd_D_tilde)]
    #v_tilde <- res_svd$v[, 1:length(svd_D_tilde)]
    #x$B <- u_tilde %*% diag(svd_D_tilde, nrow = length(svd_D_tilde)) %*% t(v_tilde)
    
    x$h_MX <- - 1/x$N * pm(x$sqrt_M_inv, t(x$x), na.rm = x$na.rm)
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
    x$QMX <- - 1/x$N * t(x$Q) %*% t(P) 
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
    x$h_K <- -1/x$N * x$K

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
    
    x$QK <- -1/x$N * t(x$Q_B) %*% x$K
    x$gamma <- 0.9 * 2 / x$D_B[1]
    x$QKMQ <- t(x$Q_B) %*% x$KM %*% x$Q_B
    
  }
  NextMethod()  
}
