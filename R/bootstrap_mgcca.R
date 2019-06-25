bootstrap = function(n, J, A, object, W, ndim){
  ind   = sample(n,replace=TRUE)
  Aboot = lapply(A, function(x) x[ind, ])

  # check for constant variables (variance equal to 0) 
  varnull = rep(0, J)
  for (xx in 1:J) varnull[xx] = length(which(apply(Aboot[[xx]], 2, var)==0))
  
  if (sum(varnull) ==0){
    result.mgcca <- mgcca_array(A = Aboot, C = object$C, tau = object$tau, ncomp = object$ncomp, scheme = object$scheme, 
                                scale = object$scale, verbose = FALSE, init = object$init, method_scale = object$method_scale)
    ########################################################################
    # Test on correlation                                                  #
    # Construction of the RGCCA outer weight vectors from boostrap sample  #
    ########################################################################
    Astar = lapply(result.rgcca$a, function(x) x[, ndim])
    Ystar = matrix(0, n, J)
    for (k in 1:J) {
      x = 1
      # Test on the sign of the correlation
      if (length(W[[k]])>1) {if (cor(W[[k]],Astar[[k]]) < 0){x  = -1}}
      else if (length(W[[k]])==1) {if (W[[k]]*Astar[[k]] < 0){x = -1}}
      Astar[[k]] = x*Astar[[k]]
      Ystar[, k] = x*result.rgcca$Y[[k]][, ndim]
    }
    MAT_COR = as.vector(cor(Ystar))
  }
  return(list(MAT_COR = MAT_COR, Astar = Astar))
}
  