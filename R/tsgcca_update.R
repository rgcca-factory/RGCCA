tsgcca_update = function(A, A_m, a, factors, weights, Y, g, dg, C,
                        ranks = rep(1, length(A)), bias = T,
                        sparsity = rep(1, length(A))) {
  ### Get useful constants
  J      <- length(A)
  DIM    <- lapply(A, dim)
  LEN    <- unlist(lapply(DIM, length))
  B_2D   <- which(LEN <= 2)   # Store which blocks are 2D
  n      <- DIM[[1]][1]
  Z      <- matrix(0, n, J)
  pjs    <- sapply(DIM, function(x) prod(x[-1]))

  ### Compute update
  for (j in 1:J) {
    const <- sparsity * sqrt(pjs)

    # Apply the derivative on the current variables
    dgx    = dg(cov2(Y[, j], Y, bias = bias))
    dgx    = matrix(rep(dgx, n), n, J, byrow = TRUE)
    Z[, j] = rowSums(matrix(rep(C[j, ], n), n, J, byrow = TRUE) * dgx * Y)

    if (j %in% B_2D) {
      # Deal with matrix blocks
      a[[j]] = t(A[[j]]) %*% Z[, j]
      a[[j]] = soft.threshold(a[[j]], const[j])
      a[[j]] = a[[j]] / norm(a[[j]], type = "2")
      Y[, j] = A[[j]] %*% a[[j]]

    } else {
      # Deal with tensor blocks
      for (d in 1:(LEN[j] - 1)) {
        other_factors = kron_prod(factors[[j]], d)
        R             = drop(weights[[j]] * sapply(1:ranks[j], function(r) prod(
          sapply(factors[[j]][-d], function(x) sum(abs(x[, r])))
        )))
        factors[[j]][[d]] = matrix(t(other_factors) %*% t(A_m[[j]]) %*% Z[, j], ncol = ranks[j])
        factors[[j]][[d]] = soft.threshold2(factors[[j]][[d]], const[j], diag(R, nrow = ranks[j]))
        factors[[j]][[d]] = apply(factors[[j]][[d]], 2, function(x) {
          mu = norm(x, type = "2")
          if (mu == 0) return(x)
          x / mu
        })

        a[[j]]            = weighted_kron_sum(factors[[j]], weights[[j]])
        Y[, j]            = A_m[[j]] %*% a[[j]]

        dgx              = dg(cov2(Y[, j], Y, bias = bias))
        dgx              = matrix(rep(dgx, n), n, J, byrow = TRUE)
        Z[, j]           = rowSums(
          matrix(rep(C[j, ], n), n, J, byrow = TRUE) * dgx * Y)
      }
      # Update weights
      if (ranks[j] > 1) {
        other_factors = list_khatri_rao(factors[[j]])
        Q             = t(other_factors) %*% t(A_m[[j]]) %*% Z[, j]
        R             = drop(sapply(1:ranks[j], function(r) sum(abs(other_factors[, r]))))
        Q             = soft.threshold3(Q, const[j], R)
        weights[[j]]  = Q / (sqrt(ranks[j]) * norm(Q, type = "2"))

        a[[j]]            = weighted_kron_sum(factors[[j]], weights[[j]])
        Y[, j]            = A_m[[j]] %*% a[[j]]
      }
    }
  }

  return(list(factors = factors, weights = weights, a = a, Y = Y))
}
