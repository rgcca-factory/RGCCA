update_mgcca = function(A, A_m, a, factors, XtX, Y, g, dg, C,
                        ranks = rep(1, length(A)), bias = T) {

  list_khatri_rao <- function(factors) {
    Reduce("khatri_rao", rev(factors))
  }

  kron_sum <- function(factors) {
    apply(list_khatri_rao(factors), 1, sum)
  }

  inv_sqrtm = function(M){
    eig        = eigen(M)
    M_inv_sqrt = eig$vectors %*% diag(eig$values^(-1/2)) %*% t(eig$vectors)
    return(M_inv_sqrt)
  }

  kron_sum_mq <- function(factors, q) {
    apply(list_khatri_rao(lapply(factors, function(x) {
      x[, -q]
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
  n   = DIM[[1]][1]
  Z   = matrix(0, n, J)

  ### Compute update
  for (j in 1:J) {
    # Apply the derivative on the current variables
    dgx    = dg(cov2(Y[, j], Y, bias = bias))
    dgx    = matrix(rep(dgx, n), n, J, byrow = TRUE)
    Z[, j] = rowSums(matrix(rep(C[j, ], n), n, J, byrow = TRUE) * dgx * Y)

    if (j %in% B_2D) {
      Az = t(A[[j]]) %*% Z[, j]
      M  = ginv(XtX[[j]])
      a[[j]] = drop(1/sqrt(t(Az) %*% M %*% Az)) * M %*% Az
      Y[, j] = A[[j]] %*% a[[j]]

    } else {
      for (r in 1:ranks[j]) {
        for (d in 1:(LEN[j] - 1)) {
          other_factors = kron_prod_q(factors[[j]], mode = d, q = r)
          Az   = t(A_m[[j]]) %*% Z[, j]
          Mqmq = t(other_factors) %*% XtX[[j]]
          Mq   = inv_sqrtm(Mqmq %*% other_factors)
          Wmq  = list_khatri_rao(lapply(factors[[j]], function(x) x[, -r, drop = F]))
          wmq  = apply(Wmq, 1, sum)
          x0   = norm(wmq, type = "2")
          wmq  = wmq / x0
          cmq  = t(wmq) %*% XtX[[j]] %*% wmq
          tmp  = Mq %*% Mqmq %*% Wmq
          pi   = tmp %*% ginv(crossprod(tmp)) %*% t(tmp)
          tmp  = Mq %*% t(other_factors) %*% Az
          lambda = sqrt(drop(t(tmp) %*% (diag(nrow(pi)) - pi) %*% tmp + (t(Az) %*% wmq)^2 / cmq))
          factors[[j]][[d]][, r] = Mq %*% (diag(nrow(pi)) - pi) %*% tmp / lambda
          x    = drop(t(Az) %*% wmq / (lambda * cmq))
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
