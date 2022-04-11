#' Function to perform the update of SGCCA variables
#'
#' @noRd
sgcca_update <- function(A, a, Y, bias, na.rm, const, J, n, dg, C) {
  Z <- matrix(0, n, J)

  for (b in seq_len(J)) {
    dgx <- dg(cov2(Y[, b], Y, bias = bias))
    CbyCovq <- drop(C[b, ] * dgx)
    Z[, b] <- Y %*% CbyCovq
    a[[b]] <- pm(t(A[[b]]), Z[, b], na.rm = na.rm)
    a[[b]] <- soft_threshold(a[[b]], const[b])
    Y[, b] <- pm(A[[b]], a[[b]], na.rm = na.rm)
  }

  return(list(a = a, Y = Y))
}
