define_M_regularisation <- function(M_regularisation, n_way, tau, A, A_m, n, p = NULL, K = NULL, J = NULL, M_K = NULL, M_J = NULL, Proj_K = NULL, Proj_J = NULL) {


  if (n_way){
    switch(M_regularisation,
     ###############################
     ##    non_kronecker_RGCCA    ##
     ###############################
     "non_kronecker_RGCCA" = 
     {
        if (tau == 1){
          P            = (n^(-1/2)) * A_m
          return(list(P = P, M_inv = NULL, tau_l = 1))
        }else{
          #Tau optimal or already define
          if (!is.numeric(tau)){
            tau_l = apply(A, 3, tau.estimate)
          }else{
            tau_l = rep(tau, K)
          }
          M_J_sqrt_inv = sapply(1:K, function(x) sqrtMatrice(tau_l[x] * diag(J) + ((1 - tau_l[x])/(n)) * (t(A[, , x]) %*% A[, , x]))$Minv_sqrt, simplify = "array")
          P            = (n^(-1/2)) * sapply(1:K, function(x) A[, , x] %*% M_J_sqrt_inv[, , x], simplify = "array")
          P            = t(apply(P, 1, c))
        }
        return(list(P = P, M_inv = NULL, tau_l = tau_l, M_J_sqrt_inv = M_J_sqrt_inv))
      },
     ###############################
     ##      kronecker_RGCCA      ##
     ###############################
     "kronecker_Identity_RGCCA" = 
     {
        P = (n^(-1/2)) * A_m
        return(list(P = P, M_inv = NULL, tau_l = 1))
      },
     ###############################
     ##      kronecker_RGCCA      ##
     ###############################
     "kronecker_specification" = 
     {
        if (!is.null(Proj_J) && !is.null(Proj_K)){
          M_J_sqrt_inv = sqrtMatrice(t(Proj_J) %*% M_J %*% Proj_J)$Minv_sqrt
          M_K_sqrt_inv = sqrtMatrice(t(Proj_K) %*% M_K %*% Proj_K)$Minv_sqrt
        }else{
          M_J_sqrt_inv = sqrtMatrice(M_J)$Minv_sqrt
          M_K_sqrt_inv = sqrtMatrice(M_K)$Minv_sqrt
        }
        P            = sapply(1:K, function(x) A[, , x] %*% M_J_sqrt_inv, simplify = "array")
        P            = aperm(a = sapply(1:J, function(x) P[, x, ] %*% M_K_sqrt_inv, simplify = "array"), perm = c(1, 3, 2))
        P            = (n^(-1/2)) * t(apply(P, 1, c))
        return(list(P = P, M_inv = NULL, tau_l = NULL, M_J_sqrt_inv = M_J_sqrt_inv, M_K_sqrt_inv = M_K_sqrt_inv))
      }
    )
  }else{
    if(tau != 1){
      if (!is.numeric(tau)) tau = tau.estimate(A)
      M_inv = ginv(tau * diag(p) + ((1 - tau)/(n)) * (t(A) %*% A))
    }else{
      M_inv = NULL
    }
    return(list(P = NULL, M_inv = M_inv, tau_l = tau))
  }
}