deflation = function(X, y, deflation_mode = NULL, bb = NULL, cc = NULL, ndefl = NULL){
  # Computation of the residual matrix R
  # Computation of the vector p.
  Proj_J = Proj_K = p = NULL
  if (length(dim(X)) > 2){
    if(is.null(deflation_mode)){
      p = apply(X, c(2, 3), miscrossprod, y)/as.vector(crossprod(y))
      R = X - aperm(a = sapply(y, function(x) x * p, simplify = "array"), perm = c(3, 1, 2))
      p = c(p)
    }else if (deflation_mode == "dim_2"){
      J      = dim(X)[2]
      K      = dim(X)[3]
      n_vect = J - ndefl
      Proj_J = lapply(1:ndefl, function(x) (diag(J) - bb[, x] %*% t(bb[, x]) / drop(t(bb[, x]) %*% bb[, x])) )
      Proj_J = Reduce("%*%", Proj_J)
      Proj_J = svd(x = Proj_J, nu = n_vect)$u
      Proj_K = diag(K)
      R      = sapply(1:K, function(x) X[, , x] %*% Proj_J, simplify = "array")
    }else if (deflation_mode == "dim_3"){
      stop(paste("Deflation mode", deflation_mode, "is not implemented yet.", sep = " "))
    }
  }else{
    p = apply(t(X),1,miscrossprod,y)/as.vector(crossprod(y))
    R = X - y%*%t(p)
  }
  return(list(p=p,R=R,Proj_J=Proj_J,Proj_K=Proj_K))
}