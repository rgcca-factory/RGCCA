
deflation <- function(X, y, left = TRUE){
  # Computation of the residual X array (can be a tensor)
  # Computation of the vector p.
  if (left) {
    if (length(dim(X)) > 2){    # Matrix case could be handled the same way
      X_m = matrix(as.vector(X), nrow = dim(X)[1])
      p   = apply(t(X_m), 1, miscrossprod, y) / drop(crossprod(y))
      R_m = X_m - pm(y,t(p))
      R   = array(as.vector(R_m), dim = dim(X))
    }else{
      p = apply(t(X),1,miscrossprod,y)/as.vector(crossprod(y))
      R = X - pm(y,t(p))
    }
  } else {
    p <- pm(X, y) / as.vector(crossprod(y))
    R <- X - tcrossprod(p, y)
  }

  return(list(p=p,R=R))
}
