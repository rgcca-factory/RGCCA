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
  if (is.null(x$f) && is.null(x$e) && is.null(x$h_tilde)) {
    NextMethod()
  } else {
    if (x$algo == 1) {
      v <- pm(x$sqrt_M, x$a, na.rm = x$na.rm)
      x$a <- (x$sqrt_M_inv %*% x$f + x$a) / drop(sqrt(
        crossprod(x$f + v)))
    } else if (x$algo == 2) {
      x$a <- t(sweep(t(x$a_MQ), 1, (x$d + 2 * x$mu)**(-1), "*")) %*% x$e
    } else if (x$algo == 3) {
      x$a <- t(sweep(t(x$a_MQ), 1, (x$eigen_val_Bplus1 + 2 * x$mu - 1)**(-1), "*")) %*% x$h_tilde 
    }
    
    x$Y <- pm(x$x, x$a, na.rm = x$na.rm)
    return(x)
  }
}
