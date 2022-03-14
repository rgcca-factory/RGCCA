rgcca_init <- function(A, init, bias, na.rm, tau, pjs, which.primal,
                       which.dual, J, n) {
  a <- alpha <- M <- Minv <- K <- list()
  Y <- matrix(0, n, J)

  if (init == "svd") {
    for (j in which.primal) {
      a[[j]] <- initsvd(A[[j]])
    }
    for (j in which.dual) {
      alpha[[j]] <- initsvd(A[[j]])
      K[[j]] <- pm(A[[j]], t(A[[j]]), na.rm = na.rm)
    }
  } else if (init == "random") {
    for (j in which.primal) {
      a[[j]] <- rnorm(pjs[j]) # random initialisation
    }

    for (j in which.dual) {
      alpha[[j]] <- rnorm(n)
      K[[j]] <- pm(A[[j]], t(A[[j]]), na.rm = na.rm)
    }
  }

  N <- ifelse(bias, n, n - 1)
  for (j in which.primal) {
    ifelse(tau[j] == 1,
           yes = {
             a[[j]] <- drop(1 / sqrt(t(a[[j]]) %*% a[[j]])) * a[[j]]
             Y[, j] <- pm(A[[j]], a[[j]], na.rm = na.rm)
           },
           no = {
             M[[j]] <- ginv(tau[j] * diag(pjs[j]) + ((1 - tau[j])) * 1 / N *
                              (pm(t(A[[j]]), A[[j]], na.rm = na.rm)))
             a[[j]] <- drop(1 / sqrt(t(a[[j]]) %*% M[[j]] %*% a[[j]])) *
               (M[[j]] %*% a[[j]])
             Y[, j] <- pm(A[[j]], a[[j]], na.rm = na.rm)
           }
    )
  }
  for (j in which.dual) {
    ifelse(tau[j] == 1,
           yes = {
             alpha[[j]] <- drop(1 / sqrt(t(alpha[[j]]) %*% K[[j]] %*%
                                           alpha[[j]])) * alpha[[j]]
             a[[j]] <- pm(t(A[[j]]), alpha[[j]], na.rm = na.rm)
             Y[, j] <- pm(A[[j]], a[[j]], na.rm = na.rm)
           },
           no = {
             M[[j]] <- tau[j] * diag(n) + ((1 - tau[j])) * 1 / N * K[[j]]
             Minv[[j]] <- ginv(M[[j]])
             alpha[[j]] <- drop(1 / sqrt(t(alpha[[j]]) %*%
               M[[j]] %*% K[[j]] %*% alpha[[j]])) * alpha[[j]]
             a[[j]] <- pm(t(A[[j]]), alpha[[j]], na.rm = na.rm)
             Y[, j] <- pm(A[[j]], a[[j]], na.rm = na.rm)
           }
    )
  }
  return(list(a = a, alpha = alpha, Y = Y, M = M, Minv = Minv, K = K))
}
