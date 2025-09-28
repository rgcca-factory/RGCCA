#' The function initsvd() is called by rgccad() and does not have to be used by
#' the user. initsvd() initializes block weight vectors based on Singular Value
#' Decomposition (SVD). Missing values are imputed by colmeans.
#' @param X  A matrix with n rows and p columns
#' @param dual A logical value. dual = TRUE enables a dual initialization (i.e.
#' the first left singular vector is used if n<p and the first right singular
#' vector is used otherwise.
#' @param rank An integer giving the number of columns extracted using the SVD.
#' @return A vector of initialization
#' @title Initialization of the T/S/RGCCA algorithm by Singular Value
#' Decomposition
#' @noRd
initsvd <- function(X, dual = TRUE, rank = 1) {
  if (any(is.na(X))) {
    indNA <- which(is.na(X), arr.ind = TRUE)
    vecMeans <- colMeans(X, na.rm = TRUE)
    X[indNA] <- vecMeans[indNA[, 2]]
  }

  n <- NROW(X)
  p <- NCOL(X)

  if (dual) {
    ifelse(n >= p,
      return(svd(X, nu = 0, nv = rank)$v),
      return(svd(X, nu = rank, nv = 0)$u)
    )
  } else {
    return(svd(X, nu = 0, nv = rank)$v)
  }
}
