mgccak_array <- function (A, A_m = NULL, C, tau = "optimal", scheme = "centroid", center = FALSE, scale = FALSE, 
                    verbose = FALSE, init="svd", bias = TRUE, tol = 1e-8, M_J = NULL, M_K = NULL, Proj_J = NULL, Proj_K = NULL,
                    M_regularisation) {
  

  # initialization
  L      <- length(A)

  #List of 3D Tensors and 2D matrix
  DIM    <- lapply(A, dim)
  LEN    <- unlist(lapply(DIM, length))
  B_3D   <- which(LEN == 3)
  B_2D   <- which(LEN == 2)
  B_0D   <- which(LEN == 0)
  
  #Convert vectors to one-column matrices
  if (length(B_0D) !=0){
    for (i in B_0D){
      A[[i]]   = as.matrix(A[[i]])
      DIM[[i]] = dim(A[[i]])
    }
    B_2D = c(B_2D, B_0D)
  }

  #Dimensions of each bloc
  Jls <- unlist(lapply(DIM, function(x) x[2]))
  Kls <- unlist(lapply(DIM, function(x) x[3]))
  pjs <- unlist(lapply(1:L, function(x) ifelse(x %in% B_3D, {return(Jls[x]*Kls[x])}, {return(Jls[[x]])})))
  n   <- DIM[[1]][1]
  Y   <- matrix(0, n, L)
  
  #Scaling
  if (scale == TRUE || center == TRUE) {
    A   = lapply(1:L, function(x) ifelse((x %in% B_3D), {return(scale3_array(A[[x]], bias = bias))}, {return(scale2(A[[x]], bias = bias)/sqrt(NCOL(x)))}))
  }
  
  #Matricization
  if(is.null(A_m)){
    A_m = lapply(1:L, function(x) ifelse(x %in% B_3D, {return(t(apply(A[[x]], 1, c)))}, {return(A[[x]])}))
  }
  
  a <- b <- c <- M_reg <- M_inv <- P <- tau_l <- list()
    
  ####################################
  # Initialisation and normalisation #
  ####################################
  
  #Initialization of vector a
  if (is.list(init)) {
    if ( !(is.null(Proj_J) | is.null(Proj_K)) ){
      a = mapply(function(x, y, z) list(t(t(x) %*% (y %x% z))), init, Proj_K, Proj_J)
    }else{
      a = init
    }
  } else if (init=="svd") {
    #SVD Initialisation of a_j
    for (j in 1:L){
      a[[j]] <- svd(A_m[[j]],nu=0,nv=1)$v
    }
  } else if (init=="random") {
    #Random Initialisation of a_j
    for (j in 1:L){
      a[[j]] <- rnorm(pjs[j])
    }
  } else {
    stop("init should be either random or by SVD.")
  }
  #Initialization of vector Y
  for (j in 1:L) Y[, j] <- (n^(-1/2)) * A_m[[j]] %*% a[[j]]
  
  N = ifelse(bias, n, n-1)  

  #Determination of the M regularization matrix
  for (j in 1:L){
    M_reg[[j]] = define_M_regularisation(M_regularisation = M_regularisation[j], 
                                         n_way            = (j %in% B_3D), 
                                         tau              = tau[j], 
                                         A                = A[[j]], 
                                         A_m              = A_m[[j]],
                                         n                = n, 
                                         p                = pjs[j], 
                                         K                = Kls[j], 
                                         J                = Jls[j], 
                                         M_K              = M_K[[j]], 
                                         M_J              = M_J[[j]],
                                         Proj_K           = Proj_K[[j]],
                                         Proj_J           = Proj_J[[j]])
    P[[j]]     = M_reg[[j]]$P
    M_inv[[j]] = M_reg[[j]]$M_inv
    tau_l[[j]] = M_reg[[j]]$tau_l
  }
  #Initialize other parameters
  iter  = 1
  crit  = numeric()
  Z     = matrix(0, NROW(A[[1]]), L)
  a_old = a

  #Compute derivate function depending on the "g" function chosen
  if (mode(scheme) == "function"){
    g  <- scheme
    dg = derivation(scheme)
  }else{
    g  <- function(x) switch(scheme,horst=x,factorial=x**2,centroid=abs(x))
    dg = NULL
  }
  #MGCCA algorithm
  repeat {
    for (j in 1:L){
      #Apply the derivate on the current variables
      dgx    = update_dgx(scheme, Y, dg, n, L, j)
      Z[, j] = rowSums(matrix(rep(C[j, ], n), n, L, byrow = TRUE) * dgx * Y)
      #Tensors
      if (j %in% B_3D){
        Q      = matrix(t(Z[, j]) %*% P[[j]], nrow = Kls[j], ncol = Jls[j], byrow = T)
        SVD    = svd(x = Q, nu = 1, nv = 1)
        c[[j]] = SVD$u
        b[[j]] = SVD$v
        Y[, j] = P[[j]] %*% (c[[j]] %x% b[[j]])
      }else{
        #Marices
        ifelse(tau[j]==1, {
          a[[j]] = drop(1/sqrt(t(Z[, j]) %*% A[[j]] %*% t(A[[j]]) %*% Z[, j])) * (t(A[[j]]) %*% Z[, j])
          Y[, j] = A[[j]] %*% a[[j]]} 
        , {
          a[[j]] = drop(1/sqrt(t(Z[, j]) %*% A[[j]] %*% M_inv[[j]] %*% t(A[[j]]) %*% Z[, j])) * (M_inv[[j]] %*% t(A[[j]]) %*% Z[, j])
          Y[, j] = A[[j]] %*% a[[j]]})
      }
    }
    #Store previous criteria
    crit[iter] <- sum(C*g(cov2(Y, bias = bias)))
    
    #After the first iteration
    if (iter > 1){
      #Change sign if needed
      for (j in B_3D){
        c[[j]] = drop(sign(cor(c[[j]], c_old[[j]])))*c[[j]]
        b[[j]] = drop(sign(cor(b[[j]], b_old[[j]])))*b[[j]]
      }
      #Compute criterion evolution
      crit_DIF          = crit[iter]-crit_old
      crit_DIF_relative = crit_DIF/crit_old
      #Compute stopping criteria
      abc               = c(Reduce("c", a[B_2D])    , Reduce("c", b[B_3D])    , Reduce("c", c[B_3D]))
      abc_old           = c(Reduce("c", a_old[B_2D]), Reduce("c", b_old[B_3D]), Reduce("c", c_old[B_3D]))
      stopping_criteria = c(sqrt(drop(crossprod(abc - abc_old))/drop(crossprod(abc_old))), crit_DIF_relative)
      #Display state of algorithm
      if (verbose){
        cat(" Iter: "        , formatC(iter,width=3, format="d"),
            " Fit: "         , formatC(crit[iter], digits=8, width=10, format="e"),
            " Dif: "         , formatC(crit_DIF, digits=8, width=10, format="e"),
            " Relative Dif: ", formatC(crit_DIF_relative, digits=8, width=10, format="e"),
            " abc Relative Dif: ", formatC(stopping_criteria[1], digits=8, width=10, format="e"),
            "\n")
      }
      #Criterion must increase
      if ( crit_DIF < -tol) stop("convergence error!!")
      #Stop if all criterion are satisfied
      if ( any(stopping_criteria < tol) | (iter > 1000)) break
    }    
    #Store previous variables
    crit_old <- crit[iter]
    a_old    <- a
    b_old    <- b
    c_old    <- c
    iter     <- iter + 1
  }

  #Inverse change of variables if needed
  for (j in B_3D){
    switch(M_regularisation[j],
     ###############################
     ##    non_kronecker_RGCCA    ##
     ###############################
     "non_kronecker_RGCCA" = 
     {
        if (tau[j] == 1){
          a[[j]] = c[[j]] %x% b[[j]]
        }else{
          a[[j]] = lapply(1:Kls[j], function(x) drop(c[[j]][x, ])*M_reg[[j]]$M_J_sqrt_inv[, , x] %*% b[[j]])
          a[[j]] = as.matrix(Reduce("c", a[[j]]))
          Y[, j] = A_m[[j]] %*% a[[j]]
        }
      },
     ###############################
     ##      kronecker_RGCCA      ##
     ###############################
     "kronecker_Identity_RGCCA" = 
     {
        if (!is.null(Proj_J[[j]]) && !is.null(Proj_K[[j]])){
          b[[j]] = Proj_J[[j]] %*% b[[j]]
          c[[j]] = Proj_K[[j]] %*% c[[j]]
          a[[j]] = c[[j]] %x% b[[j]]
        }else{
          a[[j]] = c[[j]] %x% b[[j]]
        }
      },
     ###############################
     ##      kronecker_RGCCA      ##
     ###############################
     "kronecker_specification" = 
     {
        if (!is.null(Proj_J[[j]]) && !is.null(Proj_K[[j]])){
          b[[j]] = Proj_J[[j]] %*% M_reg[[j]]$M_J_sqrt_inv %*% b[[j]]
          c[[j]] = Proj_K[[j]] %*% M_reg[[j]]$M_K_sqrt_inv %*% c[[j]]
          a[[j]] = c[[j]] %x% b[[j]]
        }else{
          b[[j]] = M_reg[[j]]$M_J_sqrt_inv %*% b[[j]]
          c[[j]] = M_reg[[j]]$M_K_sqrt_inv %*% c[[j]]
          a[[j]] = c[[j]] %x% b[[j]]
          Y[, j] = A_m[[j]] %*% a[[j]]
        }
      }
    )       
  }
  #Final messages
  if (iter > 1000)           warning("The MGCCA algorithm did not converge after 1000 iterations.")
  if (iter < 1000 & verbose) cat("The MGCCA algorithm converged to a solution of the stationary equations after", iter-1, "iterations \n")
  if (verbose)               plot(crit[1:iter], xlab = "iteration", ylab = "criteria")

  AVEinner <- sum(C * cor(Y)^2/2)/(sum(C)/2)
  
  result <- list(Y         = Y, 
                 a         = a, 
                 b         = b, 
                 c         = c, 
                 crit      = crit,
                 AVE_inner = AVEinner, 
                 C         = C, 
                 tau       = tau_l, 
                 scheme    = scheme)
  
  return(result)
}
