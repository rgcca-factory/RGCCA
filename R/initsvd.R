# The function initsvd() is called by rgccad() and does not have to be used by
# the user. initsvd() initializes block weight vectors based on Singular Value
# Decomposirion (SVD). Missing values are imputed by colmeans.
# @param X  A matrix with n rows and p columns
# @param dual A logical value. dual = TRUE enables a dual initialization (i.e.
# the first left singular vector is used is n<p and the first right singular
# vector is used otherwise.
# @return A vector of initialization
# @title Initialisation of the S/RGCCA algorithm by Singular Value Decomposition

initsvd <- function(X, dual = TRUE) {
  if (any(is.na(X))) {
    indNA <- which(is.na(X), arr.ind = TRUE)
    vecMeans <- colMeans(X, na.rm = TRUE)
    X[indNA] <- vecMeans[indNA[, 2]]
  }

  n <- NROW(X)
  p <- NCOL(X)

  if (dual) {
    ifelse(n >= p,
      return(svd(X, nu = 0, nv = 1)$v),
      return(svd(X, nu = 1, nv = 0)$u)
    )
  } else {
    return(svd(X, nu = 0, nv = 1)$v)
  }
}
