#' @importFrom MASS ginv

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
block_init.tensor_block <- function(x, init = "svd") {
  if (init == "svd") {
    x$factors <- lapply(seq_along(dim(x$x))[-1], function(m) {
      svd(apply(x$x, m, c), nu = 0, nv = x$rank)$v
    })
  } else {
    x$factors <- lapply(seq_along(dim(x$x))[-1], function(m) {
      if (m == x$mode_orth) {
        svd(matrix(
          rnorm(dim(x$x)[m] * x$rank), dim(x$x)[m]
        ), nu = x$rank, nv = 0)$u
      } else {
        matrix(rnorm(dim(x$x)[m] * x$rank), dim(x$x)[m])
      }
    })
  }
  x$weights <- rep(1 / sqrt(x$rank), x$rank)

  return(block_project(x))
}

#' @export
block_init.regularized_tensor_block <- function(x, init = "svd") {
  # Compute the highest singular value of the regularization matrix
  p <- prod(dim(x$x)[-1])
  if (p > x$n) {
    x$M <- eigen(
      pm(matrix(x$x, x$n), t(matrix(x$x, x$n)), na.rm = x$na.rm),
      symmetric = TRUE, only.values = TRUE
    )$values[1]
  } else {
    x$M <- eigen(
      pm(t(matrix(x$x, x$n)), matrix(x$x, x$n), na.rm = x$na.rm),
      symmetric = TRUE, only.values = TRUE
    )$values[1]
  }
  x$M <- x$tau + (1 - x$tau) * x$M / x$N

  # Initialize the factors and weights using the tau = 1 strategy
  x <- NextMethod()

  # Change weights to satisfy the constraints
  x$weights <- x$weights / sqrt(x$M)
  x$a <- x$a / sqrt(x$M)
  x$Y <- x$Y / sqrt(x$M)
  return(x)
}

#' @export
block_init.separable_regularized_tensor_block <- function(x, init = "svd") {
  # Compute separable estimation of the regularization matrix
  d <- length(dim(x$x)) - 1
  x$M <- estimate_separable_covariance(x$x)
  x$M <- lapply(x$M, function(y) {
    sqrt_matrix(
      x$tau^(1 / d) * diag(nrow(y)) + (1 - x$tau^(1 / d)) * y,
      inv = TRUE
    )
  })

  # Make a change of variables
  for (m in seq_len(d)) {
    x$x <- mode_product(x$x, x$M[[m]], m = m + 1)
  }

  # Initialize the factors and weights using the tau = 1 strategy
  NextMethod()
}
