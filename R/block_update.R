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
  if (x$algo == 1) {
    x$f <- pm(x$f_left, 1/x$N * grad, na.rm = x$na.rm) - 
      pm(x$f_right, x$a, na.rm = x$na.rm)
  } else if (x$algo == 3) {
    x$h_tilde <- pm(x$h_tilde_QMX, grad, na.rm = x$na.rm)
    mu_max <- 0.5 * sqrt(sum(x$h_tilde**2)) #mu_max <- 0.5 * sum(x$h_tilde**2 / (x$D + .Machine$double.eps))
    mu_min <- max(0, sapply(seq_len(x$p), function(k) {0.5 * sqrt(sum(x$h_tilde[k:x$p]**2)) - max(x$D[k:x$p])}))
    #cat("mu_min = ", max(sapply(seq_len(x$p), function(k) {0.5 * sqrt(sum(x$h_tilde[k:x$p]**2)) - max(x$D[k:x$p])})), "\n")
    L <- function(mu) {0.5 * sum(x$h_tilde**2 / (x$D + 2 * mu)) + mu}
    #cat(optimize(f=L, interval = c(0, 100000), tol = .Machine$double.eps)$minimum < mu_max, sep = " ")
    #if (!(optimize(f=L, interval = c(0, 100000))$minimum < mu_max)) {
    #  cat(optimize(f=L, interval = c(0, 100000))$minimum - mu_max, sep = "\n")
    #}
    x$mu <- optimize(f = L, interval = c(mu_min, mu_max), tol = .Machine$double.eps)$minimum
    #cat("mu = ", x$mu, "\n")
    #cat("L(mu) = ", L(x$mu), "\n")
  } else if (x$algo == 4) {
    x$h <- pm(x$h_MX, grad, na.rm = x$na.rm)
    mu_max <- 10000 #TODO find expression
    L <- function(mu) {drop(0.5 * t(x$h) %*% solve(x$B + 2 * mu * diag(nrow = x$p)) %*% x$h + mu)}
    x$mu <- optimize(f = L, interval = c(0, mu_max), tol = .Machine$double.eps)$minimum
    #plot(x$h)
    #plot(eigen(x$B + 2 * x$mu * diag(nrow = x$p), only.values = T, symmetric = T)$values)
    cat("min d_k = ", min(eigen(x$B + 2 * x$mu * diag(nrow = x$p), only.values = T)$values), "\n")
    cat("mu = ", x$mu, "\n")
    cat("L(mu) = ", L(x$mu), "\n")
  } else if (x$algo == 5) {
    x$h_tilde <- x$QMX %*% grad
    
    z_old <- x$z
    iter <- 0
    crit <- c()
    
    repeat{
      z <- z_old - x$gamma * (x$h_tilde + x$D * z_old)
      
      if (any(z != 0) && norm(z, type = "2") > 1) {
        z <- z / norm(z, type = "2")
      }
      
      crit <- c(crit, crossprod(x$h_tilde, z) + 0.5 * sum(x$D * z**2))
      
      if (iter == 10000 || norm(z - z_old, type = "2") < 1E-8) { #TODO normalize diff?
        cat("n iter = ", iter, " ; norm of z=", norm(z, type = "2"), "\n")
        break
      }
      
      iter <- iter + 1
      z_old <- z
    }
    plot(crit, col = "gold")
    x$z <- z
  }
  
  return(block_project(x))
}

#' @export
block_update.dual_ac_block <- function(x, grad) {
  if (x$algo == 3) {
    x$h <- x$h_K %*% grad
    mu_max <- 1E10
    L <- function(mu) {
      0.5 * t(x$h) %*% ginv(x$B + 2 * mu * x$KM) %*% x$h + mu
    }
    x$mu <- optimize(f = L, interval = c(0, mu_max))$minimum
  
    #x$M_grad <- 1/x$N * x$M_n_inv %*% grad
    #mu_max <- 0.5 * (1/x$N**2) * pm(t(grad), pm(x$K_M, grad, na.rm = x$na.rm), na.rm = x$na.rm)
    # L <- function(mu) {
    #   tmp <- ginv(x$O + 2 * mu * x$M_n_inv)
    #   res <- 0.5 * 1/x$N * pm(
    #     t(grad), pm(
    #       t(x$K_M), pm(
    #         tmp, 
    #         x$M_grad, na.rm = x$na.rm), na.rm = x$na.rm), na.rm = x$na.rm
    #     ) + mu
    #   return(res)
    # }
    # x$mu <- optimize(f = L, interval = c(0, 1000))$minimum
  } else if (x$algo == 5) {
    x$h_tilde <- x$QK %*% grad
    
    z_old <- x$z
    iter <- 0
    crit <- c()
    
    repeat{
      z <- z_old - x$gamma * (x$h_tilde + x$D_B * z_old)
      
      if (any(z != 0) && drop(t(z) %*% x$QKMQ %*% z) > 1) {
        y <- t(x$Q_A) %*% z
        r <- sum(x$D_A > 0)
        mu_max <- sqrt(sum(y[1:r]**2) / x$D_A[r]) - 1/x$D_A[1]
        L <- function(mu) {mu + sum(y**2 / (mu * x$D_A + 1))}
        mu <- optimize(f = L, interval = c(0, mu_max), tol = .Machine$double.eps)$minimum
        if (abs(mu_max - mu) < 1E-4) {cat("at iter ", iter, ", mu_max is really close to mu", "\n")}
        z <- x$Q_A %*% (y / (mu * x$D_A + 1))
      }
      
      crit <- c(crit, crossprod(x$h_tilde, z) + 0.5 * sum(x$D * z**2))
      
      if (iter == 10000 || norm(z - z_old, type = "2") < 1E-8) { #TODO normalize diff?
        cat(iter, "\n")
        cat("norm of z=", norm(z, type = "2"), "\n")
        break
      }
      
      iter <- iter + 1
      z_old <- z
    }
    plot(crit, col = "gold")
    x$z <- z
    
  }
  return(block_project(x))
}