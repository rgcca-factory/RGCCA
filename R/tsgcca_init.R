tsgcca_init = function(A, A_m, C, ranks = rep(1, length(A)),
                      sparsity = rep(1, length(A)), init = "svd", bias = T,
                      scheme = "centroid", tol = 1e-8) {
  ### Get useful constants
  J      <- length(A)
  DIM    <- lapply(A, dim)
  LEN    <- unlist(lapply(DIM, length))
  B_2D   <- which(LEN <= 2)   # Store which blocks are 2D

  # Dimensions of each block
  pjs <- sapply(DIM, function(x) prod(x[-1]))
  sjs <- sapply(DIM, function(x) sum(x[-1]))

  # Run MGCCA
  fit_mgcca = mgccak(A, A_m, C, scheme = scheme, init = init, bias = bias,
                     tol = tol, ranks = ranks, regularisation_matrices = NULL)

  a = list()
  factors = fit_mgcca$factors
  weights = fit_mgcca$weights

  const <- sparsity * sqrt(pjs)

  weights = lapply(1:J, function(j) if (!(j %in% B_2D)) weights[[j]] / sqrt(ranks[j]))
  for (j in 1:J) {
    if (j %in% B_2D) {
      a[[j]] = soft.threshold(fit_mgcca$a[[j]], const[j])
      a[[j]] = a[[j]] / norm(a[[j]], type = "2")
    } else {
      for (m in 1:(LEN[j] - 1)) {
        for (r in 1:ranks[j]) {
          factors[[j]][[m]][, r] = soft.threshold(factors[[j]][[m]][, r], (const[j] / (ranks[j] * abs(weights[[j]][r])))^(DIM[[j]][m + 1] / sjs[j]))
          factors[[j]][[m]][, r] = factors[[j]][[m]][, r] / norm(factors[[j]][[m]][, r], type = "2")
        }
      }
      a[[j]] = weighted_kron_sum(factors[[j]], weights[[j]])
    }
  }
  return(list(factors = factors, weights = weights, a = a))
}
