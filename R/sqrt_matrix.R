#' Compute the square root or the inverse of the square root of a
#' symmetric matrix.
#' @param X A symmetric matrix.
#' @param tol A relative tolerance to detect zero singular values.
#' @param inv A boolean indicating if the inverse of the square root must be
#' computed.
#' @noRd
sqrt_matrix <- function(X, tol = sqrt(.Machine$double.eps), inv = FALSE) {
  eig <- eigen(X, symmetric = TRUE)
  positive <- eig$values > max(tol * eig$values[1], 0)
  d <- eig$values
  if (inv) {
    d[positive] <- 1 / sqrt(d[positive])
  } else {
    d[positive] <- sqrt(d[positive])
  }
  d[!positive] <- 0
  eig$vectors %*% diag(d, nrow = length(d)) %*% t(eig$vectors)
}
