proj_l1_l2 <- function(argu, a=1){
  #Check if constraints are already satisfied
  norm2_argu = norm(argu, type = "2")
  if ( norm2_argu < 1e-32 ) stop("Norm2 of argu is too small :", norm2_argu)
  if ( sum(abs(argu/norm2_argu)) <= a ) return(list(k=NaN, lambda = 0))
  # The desired a_k cannot be null as the constraints are not already satisfied 
  # (cf. previous check). So zero values are removed.
  uneq <- argu != 0
  L    <- sum(!uneq)
  p    <- abs(argu[uneq])
  #Check if multiple maximum
  MAX  <- max(p)
  bMAX <- p == MAX
  nMAX <- sum(bMAX)
  # If there are multiple maximum value, the sparse parameter 
  # "a" must be >= sqrt(number of max)
  if (a < sqrt(nMAX)){
    warning("Impossible to project, minimum ratio is : sqrt(nMAX) = ", sqrt(nMAX), 
            ". Hence the sparse parameter is changed for sqrt(nMAX).")
    a = sqrt(nMAX)
  } 
  # If there are multiple maximum value and a = sqrt(number of max), 
  # solution is straightforward
  if (a == sqrt(nMAX)){
    if (length(bMAX) == length(p)){
      # Case where there is as many maximum as the number of elements of the 
      # vector to project "argu". Indeed in the case, 
      # ||argu||_1/||argu||_2 = sqrt(number of max)
      lambda = 0
    }else{
      # With this choice of lambda, the soft-thresholding will set all 
      # parameters to 0 except for the elements equals to the maximum value.
      a_2    = max(p[-which(bMAX)])
      lambda = (MAX + a_2)/2
    }
    return(list(k=NaN, lambda = lambda))
  }
  # If the vector to project "argu" is composed of 2 elements only, as the 
  # sparse parameter a <= 1 the desired a_k is a_2 (the lowest element) as 
  # psi(a_2) = 1 (cf. theory) and by construction a_3 = 0 and psi(0) >= a 
  # (checked previously). Moreover, these 2 elements are necessarily different 
  # (cf. conditions above)
  if (length(p) == 2){
    a_k     = min(p)
    psi_a_k = 1
    k       = 2
    lambda  = a_k - (a*sqrt((k - psi_a_k^2)/(k - a^2))-psi_a_k)*(sum(p) - k*a_k)/(psi_a_k*(k))
    return( list(k=NaN, lambda = lambda) )
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
    # NOTE : could create  : ssq   <- function(u) sum(u**2) -> not necessary 
    # could use norm(u, type = "2")^2 -> not working when u = as.numeric(0) (mean p_low is empty)
    # When p_low is empty, sum(p_low**2) = 0, which is what is wanted.
    s_low_2 = sum(p_low**2) + nb_a_k*aksq
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
  #Compute lambda
  lambda <- a_k - (a*sqrt((k - psi_a_k^2)/(k - a^2))-psi_a_k)*(s_1 + s_low_1 - k*a_k)/(psi_a_k*(k))
  return( list(k=k, lambda = lambda) )
}
