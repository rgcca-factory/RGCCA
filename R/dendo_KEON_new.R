dendo_KEON_new = function(A){
  J = length(A)
  
  KEON_1                                      = combn(x = 1:J, m = 2, FUN = function(x) mckeon(A[x]))
  KEON_1_TRI                                  = matrix(0, J, J)
  KEON_1_TRI[upper.tri(KEON_1_TRI, diag = F)] = KEON_1
  KEON_1_SYM                                  = diag(J) + KEON_1_TRI + t(KEON_1_TRI)
  
  out        = list()
  out$merge  = matrix(0, nrow = J-1, ncol = 2)
  out$height = rep(0, J - 1)
  out$order  = rep(0, J)
  out$labels = names(A)
  out$method = "McKeon"
  attr(out, "class") = "hclust"
  
  names(A) = -(1:J)
  
  i = 1
  
  L = A
  
  # browser()
  
  while (i < J){
    A = lapply(L, scale2)
    A = lapply(A, function(x) x/sqrt(NCOL(x)))
    
    J_current = length(A)
    idx_2D    = combn(x = 1:J_current, m = 2)
    FUSION    = F
    
    while (!FUSION){
      
      KEON            = apply(X = idx_2D, 2, FUN = function(x) mckeon(A[x]))
      idx_2D_max_KEON = which.max(KEON)
      idx_melt        = idx_2D[ , idx_2D_max_KEON]
      max_KEON        = max(KEON)
      
      melt_1 = as.numeric(names(A)[idx_melt[1]])
      if (melt_1 < 0){
        max_1   = 1
      }else{
        max_1   = 1 - out$height[melt_1]
      }
      
      melt_2 = as.numeric(names(A)[idx_melt[2]])
      if (melt_2 < 0){
        max_2   = 1
      }else{
        max_2   = 1 - out$height[melt_2]
      }
      
      if (any(max_KEON > c(max_1, max_2))){
        if (dim(idx_2D)[2] != 1){
          idx_2D = as.matrix(idx_2D[ ,-idx_2D_max_KEON])
        }else{
          stop("Did not converge!!")
        }
      }else{
        FUSION = T
        
        L[[J_current+1]]      = cbind(L[[idx_melt[1]]], L[[idx_melt[2]]])
        names(L)[J_current+1] = i
        out$merge[i, 1]       = melt_1
        out$merge[i, 2]       = melt_2
        out$height[i]         = 1 - max_KEON
        L[idx_melt]           = NULL
        
        i = i + 1
      }
    }
    print(i)
  }
  
  out$order = order.optimal(dist = as.dist(KEON_1_SYM), merge = out$merge)$order
  
  return(out)
}
