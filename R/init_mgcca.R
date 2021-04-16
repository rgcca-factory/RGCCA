init_mgcca = function(A, A_m, ranks = rep(1, length(A)),
                      tau = rep(1, length(A)), init = "svd", bias = T) {

  list_khatri_rao <- function(factors) {
    Reduce("khatri_rao", rev(factors))
  }

  kron_sum <- function(factors) {
    apply(list_khatri_rao(factors), 1, sum)
  }

  kron_sum_lq <- function(factors, q) {
    apply(list_khatri_rao(lapply(factors, function(x) {
      rank = ncol(x)
      x[, -(q:rank), drop = FALSE]
    })), 1, sum)
  }

  kron_prod_q <- function(factors, mode, q) {
    D = length(factors)
    Reduce("%x%", rev(lapply(1:D, function(d) {
      if (mode == d) {
        return(diag(nrow(factors[[mode]])))
      } else {
        return(factors[[d]][, q])
      }
    })))
  }

  project_factor_q <- function(factors, mode, q, M) {
    projection_matrices = list()
    res = factors[[mode]][, q]
    for (r in 1:(q - 1)) {
      res = res - factors[[mode]][, r] * drop(
        t(factors[[mode]][, r]) %*% M %*% factors[[mode]][, q] / t(factors[[mode]][, r]) %*% M %*% factors[[mode]][, r]
      )
    }
    return(res)
  }

  weighted_factor <- function(u, d, rank) {
    if (rank == 1)
      return(u)
    return((u %*% diag(d[1:rank])) / sqrt(sum(d[1:rank] ^ 2)))
  }

  ### Get useful constants
  J      <- length(A)

  # List of 2D matrix and higher order tensors
  DIM    <- lapply(A, dim)
  LEN    <- unlist(lapply(DIM, length))
  B_2D   <- which(LEN == 2)   # Store which blocks are 2D
  B_0D   <- which(LEN == 0)   # Store which blocks are 1D (stored as 0D)

  # Convert vectors to one-column matrices
  if (length(B_0D) != 0) {
    for (i in B_0D) {
      A[[i]]   = as.matrix(A[[i]])
      DIM[[i]] = dim(A[[i]])
    }
    B_2D = c(B_2D, B_0D)
  }

  # Dimensions of each block
  pjs <- sapply(DIM, function(x) prod(x[-1]))
  n   <- DIM[[1]][1]

  XtX = lapply(1:J, function(j) {
    (1 - tau[j]) * crossprod(A_m[[j]]) / (n - 1 + bias) + tau[j] * diag(pjs[j])
  })

  ### Compute factors
  a = factors = list()
  for (j in 1:J) {
    factors[[j]] = list()
    if (j %in% B_2D) {
      a[[j]] <- initsvd(A[[j]], dual = FALSE)
      a[[j]] = a[[j]] / sqrt(drop(t(a[[j]]) %*% XtX[[j]] %*% a[[j]]))
    } else {
      # Initialize weight factors
      SVD               <- svd(apply(A[[j]], 2, c), nu = 0, nv = ranks[[j]])
      factors[[j]][[1]] <- weighted_factor(SVD$v, SVD$d, ranks[j])
      for (d in 2:(LEN[[j]] - 1)) {
        factors[[j]][[d]] <- svd(apply(A[[j]], d + 1, c), nu = 0, nv = 1)$v
        factors[[j]][[d]] <- matrix(c(factors[[j]][[d]]), nrow = nrow(factors[[j]][[d]]), ncol = ranks[j])
      }

      # Change weight factors to respect orthogonality constraints
      if (ranks[j] > 1) {
        for (r in 2:ranks[j]) {
          p   = kron_prod_q(factors[[j]], mode = 1, q = r)
          M   = t(p) %*% XtX[[j]] %*% p
          factors[[j]][[1]][, r] = project_factor_q(factors[[j]], mode = 1, q = r, M = M)
          # wmq = kron_sum_lq(factors[[j]], r)
          # y   = M %*% wmq
          # factors[[j]][[1]][, r] = (diag(DIM[[j]][[2]]) - tcrossprod(y) / drop(crossprod(y))) %*% factors[[j]][[1]][, r]
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
