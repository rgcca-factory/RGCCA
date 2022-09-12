ns_mgcca_update = function(A, A_m, a, factors, weights, XtX, XtX_sing, Y, g, dg, C,
                        ranks = rep(1, length(A)), bias = T, tau = tau, tol = 1e-8) {
  ### Get useful constants
  J      <- length(A)
  DIM    <- lapply(A, dim)
  LEN    <- unlist(lapply(DIM, length))
  B_nD   <- which(LEN >= 4)   # Store which blocks are higher order tensors
  B_3D   <- which(LEN == 3)   # Store which blocks are 3D
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
      if (j %in% B_3D) { # 3D Tensors
        Q                 = matrix(t(Z[, j]) %*% A_m[[j]], nrow = DIM[[j]][3],
                                   ncol = DIM[[j]][2], byrow = T)
        SVD               = svd(x = Q, nu = ranks[[j]], nv = ranks[[j]])
        factors[[j]][[1]] = SVD$v
        factors[[j]][[2]] = SVD$u
        a[[j]]            = weighted_kron_sum(factors[[j]], weights[[j]])
        Y[, j]            = A_m[[j]] %*% a[[j]]
        Z[, j]           = rowSums(
          matrix(rep(C[j, ], n), n, J, byrow = TRUE) * dgx * Y)
      } else if (j %in% B_nD) {
      # Deal with tensor blocks
        for (d in 1:(LEN[j] - 1)) {
          Q                 = array(t(A_m[[j]]) %*% Z[, j], dim = DIM[[j]][-1])
          Q                 = unfold(Q, mode = d)
          other_factors     = list_khatri_rao(factors[[j]][-d]) %*% diag(weights[[j]], nrow = length(weights[[j]]))
          SVD               = svd(x = Q %*% other_factors, nu = ranks[j],
                                  nv = ranks[j])
          factors[[j]][[d]] = SVD$u %*% t(SVD$v)

          a[[j]]            = weighted_kron_sum(factors[[j]], weights[[j]])
          Y[, j]            = A_m[[j]] %*% a[[j]]

          dgx              = dg(cov2(Y[, j], Y, bias = bias))
          dgx              = matrix(rep(dgx, n), n, J, byrow = TRUE)
          Z[, j]           = rowSums(
            matrix(rep(C[j, ], n), n, J, byrow = TRUE) * dgx * Y)
        }
      }
      # Update weights
      W            = list_khatri_rao(factors[[j]])
      tmp          = t(Z[, j]) %*% A_m[[j]] %*% W

      lambda1 = drop(ginv(
        ((1 - tau[j]) / (n - 1 + bias) * crossprod(A_m[[j]] %*% W) + tau[j] * diag(ranks[j]))
      ) %*% t(tmp))
      lambda2 = drop(t(tmp)) / (norm(tmp, type = "2") * sqrt(XtX_sing[[j]]))
      eps <- 0.5 * drop(tmp %*% weights[[j]] + tmp %*% lambda2)

      lambda1_norm = lambda1 / (norm(lambda1, type = "2") * sqrt(XtX_sing[[j]]))
      # If lambda1 is satisfying, keep lambda1, otherwise find a point that
      # increases the criterion and satisfies the constraints between
      # lambda1 and lambda2
      if (tmp %*% lambda1_norm >= eps) {
        weights[[j]] <- lambda1_norm
      } else {
        # Otherwise, lambda1 is orthogonal to weights[[j]] so we should keep the
        # point that maximizes the criterion
        if (1 / XtX_sing[[j]] - eps^2 / (tmp %*% t(tmp)) > tol) {
          Pz = diag(ranks[j]) - crossprod(tmp) / drop(tcrossprod(tmp))
          mu = norm(Pz %*% lambda1, type = "2") / drop(sqrt(1 / XtX_sing[[j]] - eps^2 / (tmp %*% t(tmp))))
          nu = drop((tmp %*% lambda1 - eps * mu) / drop(tmp %*% t(tmp)))
          weights[[j]] = (lambda1 - nu * t(tmp)) / mu
        } else {
          weights[[j]] = lambda2
        }
      }
      weights[[j]] <- drop(weights[[j]])

      a[[j]] = weighted_kron_sum(factors[[j]], weights[[j]])
      Y[, j] = A_m[[j]] %*% a[[j]]
    }
  }

  return(list(factors = factors, weights = weights, a = a, Y = Y))
}
