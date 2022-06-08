ns_mgcca_penalized_update = function(A, A_m, a, factors, weights, XtX, XtX_sing, Y, g, dg, C,
                        ranks = rep(1, length(A)), bias = T, penalty_coef = 0) {
  ### Get useful constants
  J      <- length(A)
  DIM    <- lapply(A, dim)
  LEN    <- unlist(lapply(DIM, length))
  B_2D   <- which(LEN <= 2)   # Store which blocks are 2D
  n      <- DIM[[1]][1]
  Z      <- matrix(0, n, J)

  ### Compute update
  for (j in 1:J) {
    # Apply the derivative on the current variables
    dgx    = dg(cov2(Y[, j], Y, bias = bias))
    dgx    = matrix(rep(dgx, n), n, J, byrow = TRUE)
    Z[, j] = rowSums(matrix(rep(C[j, ], n), n, J, byrow = TRUE) * dgx * Y)

    if (j %in% B_2D) {
      # Deal with matrix blocks
      Az     = t(A[[j]]) %*% Z[, j]
      M      = ginv(XtX[[j]])
      a[[j]] = drop(1/sqrt(t(Az) %*% M %*% Az)) * M %*% Az
      Y[, j] = A[[j]] %*% a[[j]]

    } else {
      # Deal with tensor blocks
      for (r in 1:ranks[j]) {
        Wmr = list_khatri_rao(lapply(factors[[j]], function(x) x[, -r, drop = F]))
        K = XtX[[j]] %*% tcrossprod(Wmr) %*% XtX[[j]]
        alpha = norm(K, type = "F")

        for (d in 1:(LEN[j] - 1)) {
          Q = t(A_m[[j]]) %*% Z[, j]
          w = Reduce("%x%", rev(lapply(factors[[j]], function(x) x[, r, drop = FALSE])))
          wmd = Reduce("%x%", rev(lapply(factors[[j]][-d], function(x) x[, r, drop = FALSE])))
          x = Q - 2 * penalty_coef * (K %*% w - alpha * w)
          y = unfold(array(x, dim = dim(A[[j]])[-1]), mode = d) %*% wmd

          factors[[j]][[d]][, r] = y / norm(y, type = "2")

          a[[j]]            = weighted_kron_sum(factors[[j]], weights[[j]])
          Y[, j]            = A_m[[j]] %*% a[[j]]

          dgx              = dg(cov2(Y[, j], Y, bias = bias))
          dgx              = matrix(rep(dgx, n), n, J, byrow = TRUE)
          Z[, j]           = rowSums(
            matrix(rep(C[j, ], n), n, J, byrow = TRUE) * dgx * Y)
        }
      }
      # Update weights
      tmp          = t(Z[, j]) %*% A_m[[j]] %*% list_khatri_rao(factors[[j]])
      weights[[j]] = drop(tmp) / (norm(drop(tmp), type = "2") * sqrt(XtX_sing[[j]] * ranks[j]))

      a[[j]]            = weighted_kron_sum(factors[[j]], weights[[j]])
      Y[, j]            = A_m[[j]] %*% a[[j]]
    }
  }

  penalty = c()
    for (j in 1:J) {
      if (j %in% B_2D) {
        penalty <- c(penalty, 0)
      } else {
        W <- list_khatri_rao(factors[[j]])
        tmp = (t(W) %*% XtX[[j]] %*% W)^2
        diag(tmp) = 0
        penalty = c(penalty, sum(tmp))
      }
    }
  penalty <- sum(penalty)

  return(list(factors = factors, weights = weights, a = a, Y = Y, penalty = penalty))
}
