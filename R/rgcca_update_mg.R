#' Function to perform the update of RGCCA variables
#' in the multi-group framework
#'
#' @noRd
rgcca_update_mg <- function(A, na.rm, tau, dg, C, a, Y, init_object, groups = NULL) {
  J <- length(A) # number of blocks
  #nb_ind <- vapply(A, NROW, FUN.VALUE = 1L) # number of individuals per group
  p <- NCOL(A[[1]]) # number of variables per group
  
  Z <- matrix(0, p, J)
  M <- init_object$M
  
  for (i in seq_along(A)) {
    Qi <- pm(t(A[[i]]), A[[i]], na.rm = na.rm) #utiliser crossprod? (3x plus rapide)
    dgx <- dg(t(Y[, i]) %*% Y) #ou est passe phi
    CbyCovi <- drop(C[i,] * dgx)
    if (tau[i] == 1) {
      Z[, i] <- Y %*% CbyCovi
      Qz <- pm(Qi, Z[, i], na.rm = TRUE) #pourquoi TRUE #Qi is symmetric so t(Qi) = Qi
      a[[i]] <- drop(1 / sqrt(crossprod(Qz))) * Qz
      Y[, i] <- pm(Qi, a[[i]], na.rm = na.rm)
    } else {
      Z[, i] <- Y %*% CbyCovi
      Qz <- pm(Qi, Z[, i], na.rm = TRUE) #pourquoi TRUE
      a[[i]] <- drop(1 / sqrt(t(Qz) %*% M[[i]] %*% Qz)) * (M[[i]] %*% Qz)
      Y[, i] <- pm(Qi, a[[i]], na.rm = na.rm)
    }
  }
  
  return(list(a = a, Y = Y))
}
