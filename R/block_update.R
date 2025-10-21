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
block_update.tensor_block <- function(x, grad) {
  grad <- array(
    pm(t(matrix(x$x, nrow = nrow(x$x))), grad, na.rm = x$na.rm),
    dim = dim(x$x)[-1]
  )
  other_factors <- NULL
  # Update factors
  for (m in seq_along(dim(x$x)[-1])) {
    grad_m <- matrix(
      aperm(grad, c(m, seq_along(dim(grad))[-m])), nrow = dim(x$x)[m + 1]
    )
    grad_m <- grad_m %*% khatri_rao(
      Reduce(khatri_rao, rev(x$factors[-seq_len(m)])), other_factors
    )
    if (m == x$mode_orth) {
      SVD <- svd(
        grad_m %*% diag(x$lambda, nrow = x$rank), nu = x$rank, nv = x$rank
      )
      x$factors[[m]] <- SVD$u %*% t(SVD$v)
    } else {
      x$factors[[m]] <- grad_m %*% diag(x$lambda, nrow = x$rank)
      x$factors[[m]] <- apply(
        x$factors[[m]], 2, function(y) y / norm(y, type = "2")
      )
    }

    other_factors <- khatri_rao(x$factors[[m]], other_factors)
  }
  # Update lambda
  x$lambda <- t(other_factors) %*% as.vector(grad)
  x$lambda <- drop(x$lambda) / norm(drop(x$lambda), type = "2")
  return(block_project(x))
}

#' @export
block_update.regularized_tensor_block <- function(x, grad) {
  grad <- array(
    pm(t(matrix(x$x, nrow = nrow(x$x))), grad, na.rm = x$na.rm),
    dim = dim(x$x)[-1]
  )
  other_factors <- NULL
  # Update factors
  for (m in seq_along(dim(x$x)[-1])) {
    grad_m <- matrix(
      aperm(grad, c(m, seq_along(dim(grad))[-m])), nrow = dim(x$x)[m + 1]
    )
    grad_m <- grad_m %*% khatri_rao(
      Reduce(khatri_rao, rev(x$factors[-seq_len(m)])), other_factors
    ) %*% diag(x$lambda, nrow = x$rank)
    if (m == x$mode_orth) {
      SVD <- svd(grad_m, nu = x$rank, nv = x$rank)
      x$factors[[m]] <- SVD$u %*% t(SVD$v)
    } else {
      x$factors[[m]] <- apply(grad_m, 2, function(y) y / norm(y, type = "2"))
    }

    other_factors <- khatri_rao(x$factors[[m]], other_factors)
  }
  # Update lambda
  u <- drop(t(other_factors) %*% as.vector(grad))

  w_ref <- drop(ginv(
    x$tau * diag(x$rank) + (1 - x$tau) * crossprod(
      pm(matrix(x$x, x$n), other_factors, na.rm = x$na.rm)
    ) / x$N
  ) %*% u)
  w_ref_norm <- w_ref / (norm(w_ref, type = "2") * sqrt(x$M))

  w_opt <- u / (norm(u, type = "2") * sqrt(x$M))

  eps <- 0.5 * drop(t(u) %*% (x$lambda + w_opt))

  # If w_ref is satisfying, keep w_ref, otherwise find a point that
  # increases the criterion and satisfies the constraints between
  # w_ref and w_opt
  if (all(w_ref_norm == w_opt)) {
    x$lambda <- w_ref_norm
  }
  else if (t(u) %*% w_ref_norm >= eps) {
    x$lambda <- w_ref_norm
  } else {
    if (1 / x$M - eps^2 / crossprod(u) > .Machine$double.eps) {
      Pu <- diag(x$rank) - tcrossprod(u / norm(u, type = "2"))
      mu <- norm(Pu %*% w_ref, type = "2") / drop(sqrt(
        1 / x$M - eps^2 / crossprod(u)
      ))
      x$lambda <- eps / drop(crossprod(u)) * u + drop(Pu %*% w_ref) / mu
    } else { # collinearity between u and w_ref
      x$lambda <- w_opt
    }
  }
  return(block_project(x))
}
