ns_mgcca_init = function(A, A_m, ranks = rep(1, length(A)),
                      tau = rep(1, length(A)), init = "svd", bias = T,
                      kronecker_covariance = F) {
  ### Get useful constants
  J      <- length(A)
  DIM    <- lapply(A, dim)
  LEN    <- unlist(lapply(DIM, length))
  B_2D   <- which(LEN <= 2)   # Store which blocks are 2D

  # Dimensions of each block
  pjs <- sapply(DIM, function(x) prod(x[-1]))
  n   <- DIM[[1]][1]

  if (kronecker_covariance) {
    XtX = lapply(1:J, function(j) {
      if (j %in% B_2D) {
        return((1 - tau[j]) * crossprod(A_m[[j]]) / (n - 1 + bias) + tau[j] * diag(pjs[j]))
      }
      fac = estimate_kronecker_covariance(A[[j]])
      to  = tau[j] ^ (1 / length(fac))
      fac = lapply(fac, function(x) {
        (1 - to) * x + to * diag(nrow(x))
      })
      return(Reduce("%x%", rev(fac)))
    })
  } else {
    XtX = lapply(1:J, function(j) {
      (1 - tau[j]) * crossprod(A_m[[j]]) / (n - 1 + bias) + tau[j] * diag(pjs[j])
    })
  }

  ### Compute factors
  a = factors = list()
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
      SVD               <- svd(apply(A[[j]], 2, c), nu = 0, nv = ranks[[j]])
      factors[[j]][[1]] <- weighted_factor(SVD$v, SVD$d, ranks[j])
      for (d in 2:(LEN[[j]] - 1)) {
        factors[[j]][[d]] <- svd(apply(A[[j]], d + 1, c), nu = 0, nv = ranks[[j]])$v
      }

      # Change weight factors to respect orthogonality constraints
      if (ranks[j] > 1) {
        for (r in 2:ranks[j]) {
          other_factors = kron_prod_q(factors[[j]], mode = 1, q = r)
          Mqmq          = t(other_factors) %*% XtX[[j]]
          Mq            = inv_sqrtm(Mqmq %*% other_factors)
          factors[[j]][[1]][, r] = project_factor_q(factors[[j]], mode = 1, q = r, Mq = Mq, Mqmq = Mqmq)
        }
      }

      # Reweight factors to respect norm constraint
      a[[j]] <- kron_sum(factors[[j]])
      weight = drop(t(a[[j]]) %*% XtX[[j]] %*% a[[j]])
      factors[[j]][[1]] <- factors[[j]][[1]] * weight ^ (-1/2)
      a[[j]] <- kron_sum(factors[[j]])
    }
  }
  return(list(factors = factors, a = a, XtX = XtX))
}
