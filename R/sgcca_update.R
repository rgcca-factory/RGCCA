#' Function to perform the update of SGCCA variables
#'
#' @noRd
sgcca_update <- function(A, bias, na.rm, sparsity, response, disjunction,
                         dg, C, a, Y, init_object) {
  J <- length(A)
  n <- NROW(A[[1]])
  pjs <- vapply(A, NCOL, FUN.VALUE = 1L)
  const <- sparsity * sqrt(pjs)

  Z <- matrix(0, n, J)
  N <- ifelse(bias, n, n - 1)

  for (j in seq_len(J)) {
    dgx <- dg(cov2(Y[, j], Y, bias = bias))
    CbyCovq <- drop(C[j, ] * dgx)
    Z[, j] <- Y %*% CbyCovq
    a[[j]] <- pm(t(A[[j]]), Z[, j], na.rm = na.rm)
    if (isTRUE(disjunction) && (j == response)) {
      M <- ginv(1 / N * pm(t(A[[j]]), A[[j]], na.rm = na.rm))
      a[[j]] <- M %*% a[[j]] / drop(sqrt(t(a[[j]]) %*% M %*% a[[j]]))
    } else {
      a[[j]] <- soft_threshold(a[[j]], const[j])
    }
    Y[, j] <- pm(A[[j]], a[[j]], na.rm = na.rm)
  }

  return(list(a = a, Y = Y))
}
