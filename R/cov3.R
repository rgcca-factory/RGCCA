# cov3
#
# Calculates the covariance on the available data
#
# @param x a matrix x
# @param y a matrix with same size as x
# @param bias if TRUE, the estimator of variance is SS/sqrt(n-1), if FALSE, it is SS/sqrt(n)
cov3 <- function(x, y = NULL, bias = TRUE) {
  n <- NROW(x)
  if (is.null(y)) {
    x <- as.matrix(x)
    if (bias) {
      W <- matrix(1, dim(x)[1], dim(x)[2])
      W[is.na(x)] <- 0
      N <- t(W) %*% W
      C <- ((N - 1) / N) * stats::cov(x, use = "pairwise.complete.obs")
    } else {
      C <- stats::cov(x, use = "pairwise.complete.obs")
    }
  } else {
    if (bias) {
      x <- as.matrix(x)
      y <- as.matrix(y)
      W1 <- matrix(1, dim(x)[1], dim(x)[2])
      W2 <- matrix(1, dim(y)[1], dim(y)[2])
      W1[is.na(x)] <- 0
      W2[is.na(x)] <- 0
      N <- t(W1) %*% W2
      C <- ((N - 1) / N) * stats::cov(x, y, use = "pairwise.complete.obs")
    } else {
      C <- stats::cov(x, y, use = "pairwise.complete.obs")
    }
  }
  return(C)
}
