#' Function to perform the update of SGCCA variables
#'
#' @noRd
sgcca_update <- function(A, bias, na.rm, sparsity, dg, C, a, Y, init_object) {
  J <- length(A)
  n <- NROW(A[[1]])
  pjs <- vapply(A, NCOL, FUN.VALUE = 1L)
  const <- sparsity * sqrt(pjs)

  Z <- matrix(0, n, J)

  for (j in seq_len(J)) {
    dgx <- dg(cov2(Y[, j], Y, bias = bias))
    CbyCovq <- drop(C[j, ] * dgx)
    Z[, j] <- Y %*% CbyCovq
    a[[j]] <- pm(t(A[[j]]), Z[, j], na.rm = na.rm)
    a[[j]] <- soft_threshold(a[[j]], const[j])
    Y[, j] <- pm(A[[j]], a[[j]], na.rm = na.rm)
  }

  return(list(a = a, Y = Y))
}
