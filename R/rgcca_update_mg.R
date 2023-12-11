#' Function to perform the update of RGCCA variables
#' in the multi-group framework
#'
#' @noRd
rgcca_update_mg <- function(A, bias, na.rm, tau, dg, C, a, Y, L, init_object, groups = NULL) {
  J <- length(A) # number of blocks
  #nb_ind <- vapply(A, NROW, FUN.VALUE = 1L) # number of individuals per group
  p <- NCOL(A[[1]]) # number of variables per group
  
  Z <- matrix(0, p, J)
  M <- init_object$M
  
  for (j in seq_along(A)) {
    #Qj <- pm(t(A[[j]]), A[[j]], na.rm = na.rm) 
    dgx <- dg(t(L[, j]) %*% L)
    CbyCovj <- drop(C[j,] * dgx)
    if (tau[j] == 1) {
      Z[, j] <- L %*% CbyCovj
      Qz <- pm(t(A[[j]]), pm(A[[j]], Z[, j], na.rm = TRUE), na.rm = na.rm) #Qj = t(A[[j]]) %*% A[[j]] is symmetric so t(Qj) = Qj
      a[[j]] <- drop(1 / sqrt(crossprod(Qz))) * Qz
      Y[[j]] <- pm(A[[j]], a[[j]], na.rm = na.rm)
      L[, j] <- pm(t(A[[j]]), Y[[j]], na.rm = na.rm)
    } else {
      Z[, j] <- L %*% CbyCovj
      Qz <- pm(t(A[[j]]), pm(A[[j]], Z[, j], na.rm = TRUE), na.rm = na.rm)
      a[[j]] <- drop(1 / sqrt(t(Qz) %*% M[[j]] %*% Qz)) * (M[[j]] %*% Qz)
      Y[[j]] <- pm(A[[j]], a[[j]], na.rm = na.rm)
      L[, j] <- pm(t(A[[j]]), Y[[j]], na.rm = na.rm)
    }
  }
  
  return(list(a = a, Y = Y, L = L))
}
