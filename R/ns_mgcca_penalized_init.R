ns_mgcca_penalized_init = function(A, A_m, ranks = rep(1, length(A)),
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

  XtX_sing <- lapply(seq(J), function(j) {
    if (j %in% B_2D) return(NULL)
    return(svd(XtX[[j]])$d[1])
  })

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
        # A[[j]] <- array(rnorm(n = pjs[[j]]), dim = c(1, DIM[[j]][-1]))
        factors[[j]] <- lapply(seq(LEN[[j]] - 1), function(d) {
          x <- matrix(rnorm(n = DIM[[j]][d + 1] * ranks[j]), ncol = ranks[j])
          x / norm(x, type = "2")
        })
      } else {
        factors[[j]] <- lapply(seq(LEN[[j]] - 1), function(d) {
          svd(apply(A[[j]], d + 1, c), nu = 0, nv = ranks[[j]])$v
        })
      }
      a[[j]] <- list_khatri_rao(factors[[j]])

      # Initialize weights to respect norm constraint
      weights[[j]] = rep(1 / sqrt(ranks[j] * XtX_sing[[j]]), ranks[j])
      a[[j]]       = a[[j]] %*% weights[[j]]
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

  return(list(factors = factors, weights = weights, a = a, XtX = XtX, XtX_sing = XtX_sing, penalty = penalty))
}
