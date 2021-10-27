deflation <- function(X, y, na.rm = TRUE){
  # Computation of the residual matrix R
  # Computation of the vector p.
  p <- pm(t(X), y, na.rm = na.rm) / as.vector(crossprod(y))
  R <- X - tcrossprod(y, p)
  return(list(p = p, R = R))
}
