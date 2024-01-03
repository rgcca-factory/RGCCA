#' Function to perform the update of RGCCA variables
#' in the multi-group framework
#'
#' @noRd
rgcca_update_mg <- function(A, bias, na.rm, tau, dg, C, a, Y, L, init_object, groups) {
  J <- length(A) # number of blocks
  p <- NCOL(A[[1]]) # number of variables per group
  
  Z <- matrix(0, p, J)
  M <- init_object$M
  
  for (j in seq_along(A)) {
    dgx <- dg(t(L[, j]) %*% L)
    CbyCovj <- drop(C[j,] * dgx)
    if (tau[j] == 1) {
      Z[, j] <- L %*% CbyCovj
      Az <- pm(A[[j]], Z[, j], na.rm = TRUE)
      Qz <- pm(t(A[[j]]), Az, na.rm = na.rm)
      a[[j]] <- drop(1 / sqrt(crossprod(Qz))) * Qz
      Y[[j]] <- pm(A[[j]], a[[j]], na.rm = na.rm)
      L[, j] <- pm(t(A[[j]]), Y[[j]], na.rm = na.rm)
    } else {
      Z[, j] <- L %*% CbyCovj
      Az <- pm(A[[j]], Z[, j], na.rm = TRUE)
      Qz <- pm(t(A[[j]]), Az, na.rm = na.rm)
      a[[j]] <- drop(1 / sqrt(t(Qz) %*% M[[j]] %*% Qz)) * (M[[j]] %*% Qz)
      Y[[j]] <- pm(A[[j]], a[[j]], na.rm = na.rm)
      L[, j] <- pm(t(A[[j]]), Y[[j]], na.rm = na.rm)
    }
  }
  
  return(list(a = a, Y = Y, L = L))
}
