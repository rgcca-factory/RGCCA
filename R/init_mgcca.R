init_mgcca = function(A, A_m, ranks = rep(1, length(A)),
                      tau = rep(1, length(A)), init = "svd", bias = T) {

  list_khatri_rao <- function(factors) {
    Reduce("khatri_rao", rev(factors))
  }

  kron_sum <- function(factors) {
    apply(list_khatri_rao(factors), 1, sum)
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

  inv_sqrtm = function(M){
    eig        = eigen(M)
    M_inv_sqrt = eig$vectors %*% diag(eig$values^(-1/2)) %*% t(eig$vectors)
    return(M_inv_sqrt)
  }

  project_factor_q <- function(factors, mode, q, Mq, Mqmq) {
    Wmq = list_khatri_rao(lapply(factors, function(x) x[, 1:(q - 1), drop = F]))
    tmp = Mq %*% Mqmq %*% Wmq
    pi  = tmp %*% ginv(crossprod(tmp)) %*% t(tmp)
    return(Mq %*% (diag(nrow(pi)) - pi) %*% Mq %*% factors[[mode]][, q])
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
