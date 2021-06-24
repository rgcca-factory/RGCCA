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
        Az            = t(A_m[[j]]) %*% Z[, j]
        other_factors = list_khatri_rao(factors[[j]][-d]) %x% diag(DIM[[j]][[d + 1]])
        Q             = t(other_factors) %*% Az
        Q             = matrix(Q, ncol = ranks[j])
        R             = drop(weights[[j]] * sapply(1:ranks[j], function(r) prod(
          sapply((1:(LEN[j] - 1))[-d], function(m) sum(abs(factors[[j]][[m]][, r])))
        )))
        # Q                 = Q %*% diag(1 / R, nrow = ranks[j])
        tryCatch(
          {
            Q = soft.threshold2(Q, const[j], diag(R, nrow = ranks[j]))
          },
          error = function(cond) {
            print("Something went wrong")
            print(cond)
          }
        )
        # Q                 = soft.threshold2(Q, const[j], diag(R, nrow = ranks[j]))
        factors[[j]][[d]] = apply(Q, 2, function(q) {
          mu = norm(q, type = "2")
          if (mu == 0) return(q)
          q / mu
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
        # Q             = Q / R
        # Q             = soft.threshold3(Q, const[j], R)
        # weights[[j]]  = (R * Q) / (ranks[j] * norm((R * Q), type = "2"))
        Q             = soft.threshold3(Q, const[j], R)
        weights[[j]]  = Q / (sqrt(ranks[j]) * norm(Q, type = "2"))

        a[[j]]            = weighted_kron_sum(factors[[j]], weights[[j]])
        Y[, j]            = A_m[[j]] %*% a[[j]]
      }
    }
  }

  return(list(factors = factors, weights = weights, a = a, Y = Y))
}
