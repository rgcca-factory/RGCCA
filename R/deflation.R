deflation <- function(X, y, na.rm = TRUE, left = TRUE) {
  # Computation of the residual matrix R
  # Computation of the vector p.
  y_norm <- drop(crossprod(y))
  if (y_norm == 0) {
    y_norm <- 1
  }
  if (left) {
    p <- pm(t(X), y, na.rm = na.rm) / y_norm
    R <- X - tcrossprod(y, p)
  } else {
    p <- pm(X, y, na.rm = na.rm) / y_norm
    R <- X - tcrossprod(p, y)
  }
  return(list(p = p, R = R))
}
