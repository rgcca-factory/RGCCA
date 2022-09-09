ns_mgcca_init = function(A, A_m, XtX, XtX_sing, ranks = rep(1, length(A)),
                      tau = rep(1, length(A)), init = "svd", bias = T) {
  ### Get useful constants
  J      <- length(A)
  DIM    <- lapply(A, dim)
  LEN    <- unlist(lapply(DIM, length))
  B_2D   <- which(LEN <= 2)   # Store which blocks are 2D

  # Dimensions of each block
  pjs <- sapply(DIM, function(x) prod(x[-1]))
  n   <- DIM[[1]][1]

  ### Compute factors
  a = factors = weights = list()
  for (j in 1:J) {
    factors[[j]] = list()
    if (j %in% B_2D) {
      if (init == "random") {
        a[[j]] <- matrix(rnorm(n = pjs[[j]]), ncol = 1)
      } else a[[j]] <- initsvd(A[[j]], dual = FALSE)
      a[[j]] = a[[j]] / sqrt(drop(t(a[[j]]) %*% XtX[[j]] %*% a[[j]]))
    } else {
      # Initialize weight factors
      if (init == "random") {
        A[[j]] <- array(rnorm(n = pjs[[j]]), dim = c(1, DIM[[j]][-1]))
      }
      for (d in 1:(LEN[[j]] - 1)) {
        factors[[j]][[d]] <- svd(apply(A[[j]], d + 1, c), nu = 0, nv = ranks[[j]])$v
      }

      # Initialize weights to respect norm constraint
      weights[[j]] = rep(1, ranks[j])
      W            = list_khatri_rao(factors[[j]])

      lambda = weights[[j]] / sqrt(drop(
        t(weights[[j]]) %*% ((1 - tau[j]) / (n - 1 + bias) * crossprod(A_m[[j]] %*% W) + tau[j] * diag(ranks[j])) %*% weights[[j]]
      ))

      if (drop(crossprod(lambda)) <= 1 / XtX_sing[[j]]) weights[[j]] = lambda
      else weights[[j]] = weights[[j]] / (sqrt(XtX_sing[[j]]) * norm(weights[[j]], type = "2"))
      a[[j]]       = weighted_kron_sum(factors[[j]], weights[[j]])
    }
  }
  return(list(factors = factors, weights = weights, a = a))
}
