rgcca_update <- function(A, a, alpha, Y, M, K, Minv, bias, na.rm, tau,
                         which.primal, which.dual, J, n, dg, C) {
  Z <- matrix(0, n, J)

  for (j in which.primal) {
    dgx <- dg(cov2(Y[, j], Y, bias = bias))
    CbyCovj <- drop(C[j, ] * dgx)
    if (tau[j] == 1) {
      Z[, j] <- Y %*% CbyCovj
      Az <- pm(t(A[[j]]), Z[, j], na.rm = TRUE)
      a[[j]] <- drop(1 / sqrt(crossprod(Az))) * Az
      Y[, j] <- pm(A[[j]], a[[j]], na.rm = na.rm)
    } else {
      Z[, j] <- Y %*% CbyCovj
      Az <- pm(t(A[[j]]), Z[, j], na.rm = TRUE)
      a[[j]] <- drop(1 / sqrt(t(Az) %*% M[[j]] %*% Az)) * (M[[j]] %*% Az)
      Y[, j] <- pm(A[[j]], a[[j]], na.rm = na.rm)
    }
  }

  for (j in which.dual) {
    dgx <- dg(cov2(Y[, j], Y, bias = bias))
    CbyCovj <- drop(C[j, ] * dgx)
    ifelse(tau[j] == 1,
      yes = {
        Z[, j] <- Y %*% CbyCovj
        alpha[[j]] <- drop(1 / sqrt(t(Z[, j]) %*% K[[j]] %*% Z[, j])) * Z[, j]
        a[[j]] <- pm(t(A[[j]]), alpha[[j]], na.rm = na.rm)
        Y[, j] <- pm(A[[j]], a[[j]], na.rm = na.rm)
      },
      no = {
        Z[, j] <- Y %*% CbyCovj
        alpha[[j]] <- drop(
          1 / sqrt(t(Z[, j]) %*% K[[j]] %*% Minv[[j]] %*% Z[, j])
        ) * (Minv[[j]] %*% Z[, j])

        a[[j]] <- pm(t(A[[j]]), alpha[[j]], na.rm = na.rm)
        Y[, j] <- pm(A[[j]], a[[j]], na.rm = na.rm)
      }
    )
  }

  return(list(a = a, alpha = alpha, Y = Y))
}
