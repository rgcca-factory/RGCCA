deflation <- function(X, y, na.rm = TRUE, left = TRUE) {
  # Computation of the residual matrix R
  # Computation of the vector p.
  if (left) {
    p <- pm(t(X), y, na.rm = na.rm) / as.vector(crossprod(y))
    R <- X - tcrossprod(y, p)
  } else {
    p <- pm(X, y, na.rm = na.rm) / as.vector(crossprod(y))
    R <- X - tcrossprod(p, y)
  }
  return(list(p = p, R = R))
}
