ns_mgcca_update = function(A, A_m, a, factors, XtX, Y, g, dg, C,
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
          other_factors = kron_prod_q(factors[[j]], mode = d, q = r)
          Az            = t(A_m[[j]]) %*% Z[, j]
          Mqmq          = t(other_factors) %*% XtX[[j]]
          Mq            = inv_sqrtm(Mqmq %*% other_factors)

          if (ranks[j] == 1) {
            x0  = 1
            cmq = 0
            pi  = 0
          } else {
            Wmq  = list_khatri_rao(lapply(factors[[j]], function(x) x[, -r, drop = F]))
            wmq  = apply(Wmq, 1, sum)
            x0   = norm(wmq, type = "2")
            wmq  = wmq / x0
            cmq  = drop(t(wmq) %*% XtX[[j]] %*% wmq)
            pi   = construct_projector(Mq, Mqmq, Wmq)
          }
          tmp    = Mq %*% t(other_factors) %*% Az
          if (cmq < 1e-6) {
            lambda = sqrt(drop(t(tmp) %*% (diag(nrow(Mq)) - pi) %*% tmp))
            x      = 0
          } else {
            lambda = sqrt(drop(t(tmp) %*% (diag(nrow(Mq)) - pi) %*% tmp + (t(Az) %*% wmq)^2 / cmq))
            x      = drop(t(Az) %*% wmq / (lambda * cmq))
          }
          factors[[j]][[d]][, r] = Mq %*% (diag(nrow(Mq)) - pi) %*% tmp / lambda
          factors[[j]][[d]][, -r] = factors[[j]][[d]][, -r] * (x / x0)

          a[[j]]            = kron_sum(factors[[j]])
          Y[, j]            = A_m[[j]] %*% a[[j]]

          dgx              = dg(cov2(Y[, j], Y, bias = bias))
          dgx              = matrix(rep(dgx, n), n, J, byrow = TRUE)
          Z[, j]           = rowSums(
            matrix(rep(C[j, ], n), n, J, byrow = TRUE) * dgx * Y)
        }
      }
    }
  }

  return(list(factors = factors, a = a, Y = Y))
}
