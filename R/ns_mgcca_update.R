ns_mgcca_update = function(A, A_m, a, factors, weights, XtX, Y, g, dg, C,
                        ranks = rep(1, length(A)), bias = T) {
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
      for (d in 1:(LEN[j] - 1)) {
        for (r in 1:ranks[j]) {
          Az   = weights[[j]][r] * t(A_m[[j]]) %*% Z[, j]
          Mqmq = alt_prod(XtX[[j]], factors[[j]], LEN[j], d, r, side = "left")
          Mq   = inv_sqrtm(alt_prod(Mqmq, factors[[j]], LEN[j], d, r, side = "right"))

          if (ranks[j] == 1) {
            pi  = 0
          } else {
            Wmq  = list_khatri_rao(lapply(factors[[j]], function(x) x[, -r, drop = F]))
            pi   = construct_projector(Mq, Mqmq, Wmq)
          }
          tmp = Mq %*% alt_prod(Az, factors[[j]], LEN[j], d, r, side = "left")
          lambda = sqrt(drop(t(tmp) %*% (diag(nrow(Mq)) - pi) %*% tmp))
          factors[[j]][[d]][, r] = Mq %*% (diag(nrow(Mq)) - pi) %*% tmp / lambda

          a[[j]]            = weighted_kron_sum(factors[[j]], weights[[j]])
          Y[, j]            = A_m[[j]] %*% a[[j]]

          dgx              = dg(cov2(Y[, j], Y, bias = bias))
          dgx              = matrix(rep(dgx, n), n, J, byrow = TRUE)
          Z[, j]           = rowSums(
            matrix(rep(C[j, ], n), n, J, byrow = TRUE) * dgx * Y)
        }

        # Update weights
        tmp          = t(Z[, j]) %*% A_m[[j]] %*% list_khatri_rao(factors[[j]])
        weights[[j]] = drop(tmp) / norm(drop(tmp), type = "2")

        a[[j]]            = weighted_kron_sum(factors[[j]], weights[[j]])
        Y[, j]            = A_m[[j]] %*% a[[j]]

        dgx              = dg(cov2(Y[, j], Y, bias = bias))
        dgx              = matrix(rep(dgx, n), n, J, byrow = TRUE)
        Z[, j]           = rowSums(
          matrix(rep(C[j, ], n), n, J, byrow = TRUE) * dgx * Y)
      }
    }
  }

  return(list(factors = factors, weights = weights, a = a, Y = Y))
}
