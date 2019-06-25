sqrtMatrice = function(M){
  # browser()
   eig        = eigen(M)
   M_sqrt     = eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors)
   Minv_sqrt  = eig$vectors %*% diag(eig$values^(-1/2)) %*% t(eig$vectors)
   return(list(M_sqrt = M_sqrt, Minv_sqrt = Minv_sqrt))
 }