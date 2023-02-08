#' Function to perform the update of RGCCA variables
#'
#' @noRd
rgcca_update <- function(A, bias, na.rm, tau, dg, C, a, Y, init_object) {
  J <- length(A) # number of blocks
  n <- NROW(A[[1]]) # number of individuals

  Z <- matrix(0, n, J)
  K <- init_object$K
  M <- init_object$M
  Minv <- init_object$Minv

  for (j in init_object$which.primal) {
    dgx <- dg(cov2(Y[, j], Y, bias = bias))
    CbyCovj <- drop(C[j, ] * dgx)
    if (tau[j] == 1) {
      Z[, j] <- Y %*% CbyCovj
      grad <- pm(t(A[[j]]), Z[, j], na.rm = TRUE)
      if (all(grad == 0)) {
        a[[j]] <- a[[j]] * 0
      } else {
        a[[j]] <- drop(1 / sqrt(crossprod(grad))) * grad
      }
      Y[, j] <- pm(A[[j]], a[[j]], na.rm = na.rm)
    } else {
      Z[, j] <- Y %*% CbyCovj
      grad <- pm(t(A[[j]]), Z[, j], na.rm = TRUE)
      if (all(grad == 0)) {
        a[[j]] <- a[[j]] * 0
      } else {
        a[[j]] <- drop(
          1 / sqrt(t(grad) %*% M[[j]] %*% grad)
        ) * (M[[j]] %*% grad)
      }
      Y[, j] <- pm(A[[j]], a[[j]], na.rm = na.rm)
    }
  }

  for (j in init_object$which.dual) {
    dgx <- dg(cov2(Y[, j], Y, bias = bias))
    CbyCovj <- drop(C[j, ] * dgx)
    ifelse(tau[j] == 1,
      yes = {
        Z[, j] <- Y %*% CbyCovj
        alpha <- drop(1 / sqrt(t(Z[, j]) %*% K[[j]] %*% Z[, j])) * Z[, j]
        a[[j]] <- pm(t(A[[j]]), alpha, na.rm = na.rm)
        Y[, j] <- pm(A[[j]], a[[j]], na.rm = na.rm)
      },
      no = {
        Z[, j] <- Y %*% CbyCovj
        alpha <- drop(
          1 / sqrt(t(Z[, j]) %*% K[[j]] %*% Minv[[j]] %*% Z[, j])
        ) * (Minv[[j]] %*% Z[, j])

        a[[j]] <- pm(t(A[[j]]), alpha, na.rm = na.rm)
        Y[, j] <- pm(A[[j]], a[[j]], na.rm = na.rm)
      }
    )
  }

  return(list(a = a, Y = Y))
}
