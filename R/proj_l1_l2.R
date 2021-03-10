proj_l1_l2 <- function(argu, a=1){
  #Check if constraints are already satisfied
  norm2_argu = norm2(argu)
  if ( norm2_argu < 1e-32 ) stop("Norm2 of argu is too small :", norm2_argu)
  if ( sum(abs(argu/norm2_argu)) <= a ) return(list(k=NaN, lambda = 0))
  uneq <- argu != 0
  L <- sum(!uneq)
  p <- abs(argu[uneq])
  #Check if multiple maximum
  MAX <- max(p)
  bMAX <- p == MAX
  nMAX <- sum(bMAX)
  if (a < sqrt(nMAX)){
    stop("Impossible to project, minimum ratio is : ", sqrt(nMAX))
  } else if (a == sqrt(nMAX)){
    warning("a == sqrt(nMAX)")
    argu_soft        <- rep(0, length(argu))
    argu_soft[bMAX]  <- 1/sqrt(nMAX)
    return(list(k=NaN, lambda = MAX - 1/sqrt(nMAX)))
  }
  #Initialize parameters
  s_1 <- s_2 <- nb <- 0
  while (T) {
    N       <- length(p)
    if (N==0) {
      warning("length(p) = 0")
      break
    }
    #Choose next a_k
    if (N%%2 == 0){
      p_reduced = p[-which.min(p)]
      a_k     = ccaPP::fastMedian(p_reduced)
    }else{
      a_k     = ccaPP::fastMedian(p)
    }
    #Make a partition of list p
    p_inf_ak <- p < a_k
    p_sup_ak <- p > a_k
    p_high  = p[p_inf_ak]
    p_low   = p[p_sup_ak]
    #Evaluation decreasing rank of a_k
    nb_a_k  = sum(p == a_k)
    k       = nb + sum(p_sup_ak) + nb_a_k
    #Compute value of the constraint
    aksq <- a_k^2
    s_low_1 = sum(p_low) + nb_a_k*a_k
    s_low_2 = norm2(p_low)^2 + nb_a_k*aksq
    psi_a_k = (s_1 + s_low_1 - k*a_k)/sqrt(s_2 + s_low_2 - 2*a_k*(s_1 + s_low_1) + k*aksq)
    #Choose partition depending on the constraint
    if (psi_a_k > a){
      if ( length(p_low) == 0 ) break
      p         = p_low
    }else{
      if (length(p_high) == 0){
        break
      }else{
        a_k_1     = max(p_high)
        psi_a_k_1 = (s_1 + s_low_1 - k*a_k_1)/sqrt(s_2 + s_low_2 - 2*a_k_1*(s_1 + s_low_1) + k*a_k_1^2)
        if (psi_a_k_1 > a){
          break
        }
        p   = p_high
        nb  = k
        s_1 = s_1 + s_low_1
        s_2 = s_2 + s_low_2
      }
    }
  }
  #Compute lambda and thus the soft tresholded vector to return
  lambda <- a_k - (a*sqrt((k - psi_a_k^2)/(k - a^2))-psi_a_k)*(s_1 + s_low_1 - k*a_k)/(psi_a_k*(k))
  return( list(k=k, lambda = lambda) )
}