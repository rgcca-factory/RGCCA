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
      svd(apply(x$x, m, c), nu = 0, nv = 1)$v
    })
  } else {
    x$factors <- lapply(seq_along(dim(x$x))[-1], function(m) {
      svd(matrix(rnorm(dim(x$x)[m]), dim(x$x)[m]), nu = 0, nv = 1)$v
    })
  }
  x$weights <- rep(1 / sqrt(1), 1)

  return(block_project(x))
}
