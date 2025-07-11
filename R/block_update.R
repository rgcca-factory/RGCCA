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
    # x$f <- pm(x$f_left, 1/x$N * grad, na.rm = x$na.rm) - 
    #   pm(x$f_right, x$a, na.rm = x$na.rm)
  } else if (x$algo == 3) {
    x$h_tilde <- x$h_tilde_QMX %*% grad
    mu_max <- 0.5 * sqrt(sum(x$h_tilde**2))
    mu_min <- max(0, sapply(seq_len(x$p), function(k) {0.5 * sqrt(sum(x$h_tilde[k:x$p]**2)) - max(x$D[k:x$p])}))
    L <- function(mu) {0.5 * sum(x$h_tilde**2 / (x$D + 2 * mu)) + mu}
    x$mu <- optimize(f = L, interval = c(mu_min, mu_max), tol = .Machine$double.eps)$minimum
    
    x$a <- x$a_MQ %*% (x$h_tilde / (x$D + 2 * x$mu))
    
  } else if (x$algo == 4) {
    # x$h <- pm(x$h_MX, grad, na.rm = x$na.rm)
    # mu_max <- 10000 #TODO find expression
    # L <- function(mu) {drop(0.5 * t(x$h) %*% solve(x$B + 2 * mu * diag(nrow = x$p)) %*% x$h + mu)}
    # x$mu <- optimize(f = L, interval = c(0, mu_max), tol = .Machine$double.eps)$minimum
    # 
  } else if (x$algo == 5) {
    x$h_tilde <- x$QMX %*% grad
    
    z_old <- x$z
    iter <- 1
    crit <- c(crossprod(x$h_tilde, z_old) + 0.5 * sum(x$D * z_old**2))
    
    repeat{
      z <- z_old - x$gamma * (x$h_tilde + x$D * z_old)
      
      if (any(z != 0) && norm(z, type = "2") > 1) {
        z <- z / norm(z, type = "2")
      }
      
      crit <- c(crit, crossprod(x$h_tilde, z) + 0.5 * sum(x$D * z**2))
      
      stopping_criteria <- c(
        norm(z - z_old, type = "2"), abs(crit[iter + 1] - crit[iter]) #/ (1 + abs(crit[iter]))
      )
      
      if (any(stopping_criteria < 1e-8) || iter == 10000) {
        #cat("n iter = ", iter, " ; norm of z=", norm(z, type = "2"), "\n")
        break
      }
      
      iter <- iter + 1
      z_old <- z
    }
    #plot(crit, col = "gold")
    x$z <- z
  }
  
  x$a <- x$MQ %*% x$z
  
  return(block_project(x))
}

#' @export
block_update.dual_ac_block <- function(x, grad) {
  if (x$algo == 3) {
    x$h <- x$h_K %*% grad
    if (!is.null(x$mu)) {
      L_mu <- -0.5 * drop(t(x$h) %*% ginv(x$B + 2 * x$mu * x$KM, tol = .Machine$double.eps * x$n) %*% x$h) - x$mu
      cat("L_(mu) =", L_mu, "\n")
    } else {
      L_mu <- NULL
    }
    
    mu_max <- 1E10
    mu_max_test <- 0.5 * drop(t(x$h) %*% ginv(x$B) %*% x$h)
    
    L <- function(mu) {
      0.5 * drop(t(x$h) %*% ginv(x$B + 2 * mu * x$KM, tol = .Machine$double.eps * x$n) %*% x$h) + mu
    }
    x$mu <- optimize(f = L, interval = c(0, mu_max))$minimum
    #cat("mu = ", x$mu, " ; mu_max_test = ", mu_max_test, ifelse(x$mu > mu_max_test, yes = "NOOOOO", no = ""), "\n")
    
    L_mu_opt <- -0.5 * drop(t(x$h) %*% ginv(x$B + 2 * x$mu * x$KM, tol = .Machine$double.eps * x$n) %*% x$h) - x$mu
    if (!is.null(L_mu)) {cat("L_(mu^) =", L_mu_opt, ifelse(L_mu - L_mu_opt > 1e-8, " NOOOO: L_(mu) > L_(mu^) ", ""), "\n")}
    
    x$alpha <- - ginv(x$B + 2 * x$mu * x$KM, tol = .Machine$double.eps * x$n) %*% x$h

  } else if (x$algo == 5) {
    x$h_tilde <- x$QK %*% grad
    
    z_old <- x$z
    iter <- 0
    crit <- c(crossprod(x$h_tilde, z_old) + 0.5 * sum(x$D * z_old**2))
    
    repeat{
      z <- z_old - x$gamma * (x$h_tilde + x$D_B * z_old)
      
      if (any(z != 0) && drop(t(z) %*% x$QKMQ %*% z) > 1) {
        y <- t(x$Q_A) %*% z
        r <- sum(x$D_A > 0)
        mu_max <- sqrt(sum(y[1:r]**2) / x$D_A[r]) - 1/x$D_A[1]
        L <- function(mu) {mu + sum(y**2 / (mu * x$D_A + 1))}
        mu <- optimize(f = L, interval = c(0, mu_max), tol = .Machine$double.eps)$minimum
        if (mu_max - mu < 1E-4) {cat("at iter ", iter, ", mu_max is really close to mu", "\n")}
        z <- x$Q_A %*% (y / (mu * x$D_A + 1))
      }
      
      crit <- c(crit, crossprod(x$h_tilde, z) + 0.5 * sum(x$D * z**2))
      
      stopping_criteria <- c(
        norm(z - z_old, type = "2") < 1e-4, abs(crit[iter + 1] - crit[iter]) < 1e-6 #/ (1 + abs(crit[iter]))
      )
      
      if (all(stopping_criteria) || iter == 10000) {
        #cat("n iter = ", iter, " ; norm of z=", norm(z, type = "2"), "\n")
        break
      }

      iter <- iter + 1
      z_old <- z
    }
    #plot(crit, col = "gold")
    x$z <- z
    
    x$alpha <- x$Q_B %*% x$z
  }
  return(block_project(x))
}