#' @importFrom MASS ginv

block_project <- function(x) {
  UseMethod("block_project")
}

#' @export
block_project.block <- function(x) {
  if (any(x$a != 0)) {
    x$a <- x$a / norm(x$a, type = "2")
  }

  x$Y <- pm(x$x, x$a, na.rm = x$na.rm)
  return(x)
}

#' @export
block_project.dual_block <- function(x) {
  if (any(x$alpha != 0)) {
    x$alpha <- x$alpha / drop(sqrt(t(x$alpha) %*% x$K %*% x$alpha))
  }
  x$a <- pm(t(x$x), x$alpha, na.rm = x$na.rm)

  x$Y <- pm(x$x, x$a, na.rm = x$na.rm)
  return(x)
}

#' @export
block_project.primal_regularized_block <- function(x) {
  if (any(x$a != 0)) {
    x$a <- x$M %*% x$a / drop(sqrt(t(x$a) %*% x$M %*% x$a))
  }

  x$Y <- pm(x$x, x$a, na.rm = x$na.rm)
  return(x)
}

#' @export
block_project.dual_regularized_block <- function(x) {
  if (any(x$alpha != 0)) {
    x$alpha <- x$M %*% x$alpha / drop(sqrt(
      t(x$alpha) %*% x$M %*% x$K %*% x$alpha
    ))
  }
  x$a <- pm(t(x$x), x$alpha, na.rm = x$na.rm)

  x$Y <- pm(x$x, x$a, na.rm = x$na.rm)
  return(x)
}

#' @export
block_project.sparse_block <- function(x) {
  if (any(x$a != 0)) {
    x$a <- soft_threshold(x$a, x$const)
  }
  x$Y <- pm(x$x, x$a, na.rm = x$na.rm)
  return(x)
}

#' @export
block_project.ac_block <- function(x) {
  if (is.null(x$f) && is.null(x$h_tilde) && is.null(x$h)) { #init
    #NextMethod()
    if (any(x$a != 0) && drop(sqrt(t(x$a) %*% x$M %*% x$a)) > 1) {
      x$a <- x$a / drop(sqrt(t(x$a) %*% x$M %*% x$a))
    }
    
    if (x$algo == 5) {
      x$z <- x$a
      if (any(x$z != 0) && norm(x$z, type = "2") > 1) {
        x$z <- x$z / norm(x$z, type = "2")
      }
      # if (any(x$a != 0) && norm(crossprod(x$Q, x$sqrt_M %*% x$a), type = "2") > 1) {
      #   x$z <- crossprod(x$Q, x$sqrt_M %*% x$a) / norm(crossprod(x$Q, x$sqrt_M %*% x$a), type = "2")
      # } else {
      #   x$z <- crossprod(x$Q, x$sqrt_M %*% x$a)
      # }
    }
    
    x$Y <- pm(x$x, x$a, na.rm = x$na.rm)
    return(x)
  } else { #update
    if (x$algo == 1) {
      # v <- pm(x$sqrt_M, x$a, na.rm = x$na.rm)
      # x$a <- (x$sqrt_M_inv %*% x$f + x$a) #/ drop(sqrt(crossprod(x$f + v)))
    } else if (x$algo == 3) {
      #x$a <- t(sweep(t(x$a_MQ), 1, x$D + 2 * x$mu, "/")) %*% x$h_tilde #not faster than the line below
      #if (norm(x$sqrt_M %*% x$a, type = "2") > 1) {cat("w^T M_j w - 1 = ", t(x$a) %*% x$M %*% x$a - 1, "\n")}
    } else if (x$algo == 4) {
      # tmp <- solve(x$B + 2 * x$mu * diag(nrow = x$p))
      # #cat("norm of v = ", crossprod(tmp %*% x$h), "\n")
      # x$a <- - x$sqrt_M_inv %*% tmp %*% x$h #/ drop(sqrt(crossprod(tmp %*% x$h)))
    } else if (x$algo == 5) {
      # x$a <- x$MQ %*% x$z
      # if (norm(x$sqrt_M %*% x$a, type = "2") > 1) {cat("w^T M_j w - 1 = ", t(x$a) %*% x$M %*% x$a - 1, "\n")}
    }
    
    x$Y <- pm(x$x, x$a, na.rm = x$na.rm)
    return(x)
  }
}

#' @export
block_project.dual_ac_block <- function(x) {
  if (is.null(x$mu) && is.null(x$h_tilde)) {# init
    if (any(x$alpha != 0) && drop(t(x$alpha) %*% x$KM %*% x$alpha) > 1) {
      x$alpha <- x$alpha / sqrt(drop(t(x$alpha) %*% x$KM %*% x$alpha)) #Is this ok?
    }
    
    if (x$algo == 5) {
      x$z <- t(x$Q_B) %*% x$alpha
    } #NextMethod()
  } # else { update
    # if (x$algo == 3) {
    #   
    # 
    # } else if (x$algo == 5) {
    #   
    #   
    # }
  # }
  
  x$a <- pm(t(x$x), x$alpha, na.rm = x$na.rm)
  
  x$Y <- pm(x$x, x$a, na.rm = x$na.rm)
  return(x)
}