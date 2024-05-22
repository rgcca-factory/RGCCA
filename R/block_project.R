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
    alpha_norm <- t(x$alpha) %*% x$K %*% x$alpha
    if (alpha_norm > 0) {
      x$alpha <- x$alpha / drop(sqrt(alpha_norm))
    } else {
      x$alpha <- x$alpha * 0
    }
  }
  x$a <- pm(t(x$x), x$alpha, na.rm = x$na.rm)

  x$Y <- pm(x$x, x$a, na.rm = x$na.rm)
  return(x)
}

#' @export
block_project.primal_regularized_block <- function(x) {
  if (any(x$a != 0)) {
    a_norm <- sqrt(t(x$a) %*% x$M %*% x$a)
    if (a_norm > 0) {
      x$a <- x$M %*% x$a / drop(sqrt(a_norm))
    } else {
      x$a <- x$a * 0
    }
  }

  x$Y <- pm(x$x, x$a, na.rm = x$na.rm)
  return(x)
}

#' @export
block_project.dual_regularized_block <- function(x) {
  if (any(x$alpha != 0)) {
    alpha_norm <- t(x$alpha) %*% x$M %*% x$K %*% x$alpha
    if (alpha_norm > 0) {
      x$alpha <- x$M %*% x$alpha / drop(sqrt(alpha_norm))
    } else {
      x$alpha <- x$alpha * 0
    }
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
