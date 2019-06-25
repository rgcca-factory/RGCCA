m_bootstrap = function(n, L, A, object, W, B, C, ndim, B_3D, paired = F){
  
  if(paired){
    ind = sample(n/2, replace=TRUE)
    ind = c(ind, ind + n/2)
  }else{
    ind = sample(n,replace=TRUE)
  }
  Aboot        = rep(list(0), L)
  Aboot[B_3D]  = lapply(A[B_3D], function(x) x[ind, , ])
  Aboot[-B_3D] = lapply(A[-B_3D], function(x) as.matrix(as.matrix(x)[ind, ]))
  
  # check for constant variables (variance equal to 0) 
  varnull = rep(0, L)
  for (xx in 1:L){
    if(xx %in% B_3D){
      varnull[xx] = length(which(apply(Aboot[[xx]], c(2, 3), var)==0))
    }else{
      varnull[xx] = length(which(apply(Aboot[[xx]], 2, var)==0))
    }
  }
  
  if (sum(varnull) ==0){
    tau          = unlist(lapply(object$tau, function(x) unique(x)))
    result.mgcca = mgcca_array(A = Aboot, C = object$C, tau = tau, ncomp = ndim, 
                               scheme = object$scheme, center = object$center, scale = object$scale, bias = T, verbose = F)
    ########################################################################
    # Test on correlation                                                  #
    # Construction of the RGCCA outer weight vectors from boostrap sample  #
    ########################################################################
    Astar = lapply(result.mgcca$a, function(x) x[, ndim])
    Bstar = lapply(result.mgcca$b[B_3D], function(x) x[, ndim])
    Cstar = lapply(result.mgcca$c[B_3D], function(x) x[, ndim])
    Ystar = Reduce("cbind", result.mgcca$Y)
    
    sign_a_cor = mapply(function(x, y) ifelse(length(x) == 1, sign(x*y), sign(drop(cor(x, y)))), Astar, W)
    sign_b_cor = mapply(function(x, y) sign(drop(cor(x, y))), Bstar, B)
    sign_c_cor = mapply(function(x, y) sign(drop(cor(x, y))), Cstar, C)
    
    Astar = lapply(1:L, function(x) sign_a_cor[x]*Astar[[x]])
    Bstar = lapply(1:length(B_3D), function(x) sign_b_cor[x]*Bstar[[x]])
    Cstar = lapply(1:length(B_3D), function(x) sign_c_cor[x]*Cstar[[x]])
    
    Ystar = Ystar %*% diag(unlist(sign_a_cor))
    
    # for (k in 1:L) {
    #   if (k %in% B_3D){
    #     x = 1
    #     # Test on the sign of the correlation
    #     if (length(B[[k]])>1) {if (cor(B[[k]],Bstar[[k]]) < 0){x  = -1}}
    #     else if (length(B[[k]])==1) {if (B[[k]]*Bstar[[k]] < 0){x = -1}}
    #     Bstar[[k]] = x*Bstar[[k]]
    #     y = 1
    #     # Test on the sign of the correlation
    #     if (length(C[[k]])>1) {if (cor(C[[k]],Cstar[[k]]) < 0){y  = -1}}
    #     else if (length(C[[k]])==1) {if (C[[k]]*Cstar[[k]] < 0){y = -1}}
    #     Cstar[[k]] = y*Cstar[[k]]
    #     Ystar[, k] = y*x*result.mgcca$Y[[k]][, ndim]
    #   }else{
    #     x = 1
    #     # Test on the sign of the correlation
    #     if (length(W[[k]])>1) {if (cor(W[[k]],Astar[[k]]) < 0){x  = -1}}
    #     else if (length(W[[k]])==1) {if (W[[k]]*Astar[[k]] < 0){x = -1}}
    #     Astar[[k]] = x*Astar[[k]]
    #     Ystar[, k] = x*result.mgcca$Y[[k]][, ndim]
    #   }
    # }
    MAT_COR = as.vector(cor(Ystar))
    # Astar   = lapply(result.mgcca$a[-B_3D], function(x) x[, ndim])
  }
  return(list(MAT_COR = MAT_COR, Astar = Astar, Bstar = Bstar, Cstar = Cstar))
}
