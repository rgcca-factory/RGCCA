mgcca_array <- function(A, C = 1-diag(length(A)), tau = rep(1, length(A)), ncomp = rep(1, length(A)), scheme = "centroid", center = TRUE, scale = TRUE , 
                        init="svd", bias = TRUE, tol = 1e-8, verbose=TRUE, method_scale = "3-way", M_J = NULL, M_K = NULL, deflation_mode = NULL,
                        M_regularisation, nstart = 1, nstart_at_comp_2 = F) {
  
  if (any(ncomp < 1)) stop("Compute at least one component per block!")

  #Number of blocks
  L      = length(A)

  #List of 3D Tensors and 2D matrix
  DIM    = lapply(A, dim)
  LEN    = unlist(lapply(DIM, length))
  B_3D   = which(LEN == 3)
  B_2D   = which(LEN == 2)
  B_0D   = which(LEN == 0)
  
  #Convert vectors to one-column matrices
  if (length(B_0D) !=0){
    for (i in B_0D){
      A[[i]]   = as.matrix(A[[i]])
      DIM[[i]] = dim(A[[i]])
    }
    B_2D = c(B_2D, B_0D)
  }
  
  #Dimensions of each bloc
  Jls    = unlist(lapply(DIM, function(x) x[2]))
  Kls    = unlist(lapply(DIM, function(x) x[3]))
  pjs    = unlist(lapply(1:L, function(x) ifelse(x %in% B_3D, {return(Jls[x]*Kls[x])}, {return(Jls[[x]])})))
  nb_row = DIM[[1]][1]
  
  #Chek if M_J and M_K are well defined
  if (!is.null(M_J)){
    if (is.null(M_K)){
      stop("Please M_J is defined but not M_K")
    }else{
      for (d in B_3D){
        if(!isSquare(M_J[[d]])) stop(paste("Please M_J n ", d, "is not a square matrix"))
        if(!isSquare(M_K[[d]])) stop(paste("Please M_K n ", d, "is not a square matrix"))
        Jk = unique(dim(M_J[[d]]))
        if(Jk !=Jls[d]) stop(paste("Please M_J n ", d, "is of dimension", Jk, "x", Jk, "instead of", Jls[d], "x", Jls[d]))
        Kk = unique(dim(M_K[[d]]))
        if(Kk !=Kls[d]) stop(paste("Please M_K n ", d, "is of dimension", Kk, "x", Kk, "isntead of", Kls[d], "x", Kls[d]))
      }
    }
  }else if (!is.null(M_K)){
    stop("Please M_K is defined but not M_J")
  }

  #Check if M_regularisation is compatible with deflation_mode
  if((M_regularisation == "non_kronecker_RGCCA") && (!is.null(deflation_mode))){
    stop(paste("Deflation mode", deflation_mode, "is not compatible with regularisation", M_regularisation, sep = " "))
  }

  #Multiple starts are random
  if (nstart > 1) init = "random"
  
  #Chek if ncomp and scheme are well defined
  if (any(ncomp-pjs > 0)) stop("For each block, choose a number of components smaller than the number of variables!")
  if ( (mode(scheme) != "function") & !(scheme %in% c("horst", "factorial", "centroid")) ){
    stop("Choose one of the three following schemes: horst, centroid, factorial or design the g function")
  }

  #Display messages
  if(verbose){
    if (mode(scheme) != "function") {
      cat("Computation of the MGCCA block components based on the", scheme, "scheme \n")
    }else{
      cat("Computation of the MGCCA block components based on the g scheme \n")
    }
    if (!is.numeric(tau)) {
      cat("Optimal Shrinkage intensity paramaters are estimated \n")
    }else{
      cat("Shrinkage intensity paramaters are chosen manually \n")
    }
  }
    
  #Scaling
  if (center == TRUE || scale == TRUE){
    A   = lapply(1:L, function(x) ifelse((x %in% B_3D), 
                                         {return(scale3_array(A[[x]], bias = bias, scale = scale, method_scale = method_scale)/sqrt(NCOL(A[[x]])))}, 
                                         {return(scale2(A[[x]], bias = bias, scale = scale)/sqrt(NCOL(A[[x]])))}))
  }
  
  #Store mean and standard deviation
  A_mean        = lapply(A, function(x) attr(x, "scaled:center"))
  names(A_mean) = names(A)
  A_sd          = NULL
  if (scale){
    A_sd        = lapply(A, function(x) attr(x, "scaled:scale"))
    names(A_sd) = names(A)
  }
  
  #Matricization of tensors
  A_m = lapply(1:L, function(x) ifelse(x %in% B_3D, {return(t(apply(A[[x]], 1, c)))}, {return(A[[x]])}))

  ######################
  ### Initialization ###
  ######################
  AVE_outer        = vector()
  ndefl            = ncomp-1
  N                = max(ndefl)
  AVE_inner        = rep(NA,max(ncomp))
  R                = A
  R_m              = A_m
  data_to_deflate  = A
  
  AVE_X <- crit <- list()
  Y     <- P    <- a <- astar <- b <- c <- Proj_J <- Proj_K <- NULL

  if(is.null(deflation_mode)){blocs_to_defl = 1:L}else{blocs_to_defl = B_2D}
  if (!is.numeric(tau)) tau_mat = c()
  if (!is.matrix(tau)){
    tau = t(sapply(1:(N+1), function(x) rep(tau, 1)))
  }else{
    if (length(tau == L)) tau = t(sapply(1:(N+1), function(x) rep(tau, 1)))
  }
  
  for (d in 1:L) P[[d]]  <- a[[d]] <- astar[[d]] <- matrix(NA,pjs[[d]],N+1)
  for (d in B_3D) b[[d]] = matrix(NA, Jls[[d]], N+1)
  for (d in B_3D) c[[d]] = matrix(NA, Kls[[d]], N+1)
  for (d in 1:L)  Y[[d]] = matrix(NA, nb_row,   N+1)
  
  #####################################
  ### Lauch algorithm per component ###
  #####################################
  for (n in 1:(N+1)) {
    if (verbose) cat(paste0("Computation of the MGCCA block components #", n, " is under progress...\n"))
    nstart_comp      = ifelse(nstart_at_comp_2, 1, nstart)
    nstart_at_comp_2 = F
    if(!is.character(init)){
      cur_init = lapply(init, function(x) x[, n])
    }else{
      cur_init = init
    }
    #n_random_starts
    for (start in 1:nstart_comp){
      #MGCCA algorithm
      tmp_mgcca.result = mgccak_array(A          = R, 
                                      A_m        = R_m, 
                                      C          = C, 
                                      tau        = tau[n, ], 
                                      scheme     = scheme, 
                                      init       = cur_init, 
                                      bias       = bias, 
                                      tol        = tol, 
                                      verbose    = verbose, 
                                      M_J        = M_J, 
                                      M_K        = M_K, 
                                      Proj_J     = Proj_J, 
                                      Proj_K     = Proj_K,
                                      M_regularisation = M_regularisation)
      tmp_crit = tmp_mgcca.result$crit[length(tmp_mgcca.result$crit)]
      if(start == 1){
        mgcca.result = tmp_mgcca.result
        best_crit    = tmp_crit
      }else{
        if(tmp_crit > best_crit){
          mgcca.result = tmp_mgcca.result
          best_crit    = tmp_crit
        }
      }
    }
    #Store tau, AVE_inner, crit
    if (!is.numeric(tau)) tau_mat = rbind(tau_mat, mgcca.result$tau)
    AVE_inner[n] = mgcca.result$AVE_inner
    crit[[n]]    = mgcca.result$crit
    
    #Store Y, a, b, c
    for (d in 1:L) Y[[d]][,n] = mgcca.result$Y[ , d]
    for (d in 1:L) a[[d]][,n] = mgcca.result$a[[d]]
    for (d in B_3D) {
      b[[d]][,n] = mgcca.result$b[[d]]
      c[[d]][,n] = mgcca.result$c[[d]]
    }
    #Deflation procedure
    defla.result = defl.select(mgcca.result$Y, data_to_deflate, ndefl, n, nbloc = L, 
                                lapply(B_3D, function(x) matrix(b[[x]][, 1:n], ncol = n, byrow = F)), 
                                lapply(B_3D, function(x) matrix(c[[x]][, 1:n], ncol = n, byrow = F)), deflation_mode)
    R                              = defla.result$resdefl
    R_m                            = NULL
    if (length(blocs_to_defl) != 0) data_to_deflate[blocs_to_defl] = R[blocs_to_defl]
      
    #Store projection matrices for deflation
    for (d in blocs_to_defl) P[[d]][,n] = defla.result$pdefl[[d]]
    if (!is.null(deflation_mode)){
      for (d in B_3D){
        Proj_J[[d]] = defla.result$Proj_J[[d]]
        Proj_K[[d]] = defla.result$Proj_K[[d]]
      }
    }
    #Compute astar
    for (d in blocs_to_defl){
      if (n == 1){
        astar[[d]][,n] = mgcca.result$a[[d]]
      }else{
        astar[[d]][,n] = mgcca.result$a[[d]] - astar[[d]][,(1:n-1), drop=F] %*% drop( t(a[[d]][,n]) %*% P[[d]][,1:(n-1),drop=F] )
      }
    }
  }

  #############
  ### Names ###
  #############
  for (d in B_3D){
    rownames(b[[d]]) = dimnames(A[[d]])[[2]]
    rownames(c[[d]]) = dimnames(A[[d]])[[3]]
    rownames(Y[[d]]) = dimnames(A[[d]])[[1]]
    colnames(Y[[d]]) = paste0("comp", 1:max(ncomp))
  }
  for (d in B_2D){
    rownames(a[[d]]) = rownames(astar[[d]]) = colnames(A[[d]])
    rownames(Y[[d]]) = rownames(A[[d]])
    colnames(Y[[d]]) = paste0("comp", 1:max(ncomp))
  }

  #Average Variance Explained (AVE) per block
  for (j in 1:L) AVE_X[[j]] =  apply(cor(A_m[[j]], Y[[j]])^2, 2, mean)
  #AVE outer 
  if (N == 0) { 
    AVE_outer = sum(pjs * unlist(AVE_X))/sum(pjs)
  }else{
    outer = matrix(unlist(AVE_X), nrow = max(ncomp))
    for (j in 1:max(ncomp)) AVE_outer[j] = sum(pjs * outer[j, ])/sum(pjs)
  }

  #For each bloc define if it is Pimal or Dual
  mode = rep(0, L)
  for (j in 1:L)  mode[j] = ifelse(nb_row>=pjs[j], "Primal", "Dual") 

  ncomp_3D = ncomp
  if (length(B_2D) !=0) ncomp_3D = ncomp[-B_2D]
  
  #Remove unused components
  if (length(unique(ncomp)) != 1){
    shave.matlist <- function(mat_list, nb_cols) mapply(function(m, nbcomp) m[, 1:nbcomp, drop = FALSE], mat_list, nb_cols,SIMPLIFY=FALSE)
    shave.veclist <- function(vec_list, nb_elts) mapply(function(m, nbcomp) m[1:nbcomp], vec_list, nb_elts, SIMPLIFY=FALSE)    
    AVE_X = shave.veclist(AVE_X, ncomp)
    Y     = shave.matlist(Y, ncomp)
    a     = shave.matlist(a, ncomp)
    b     = shave.matlist(b, ncomp_3D)
    c     = shave.matlist(c, ncomp_3D) 
    astar = shave.matlist(astar, ncomp)
  }

  #tau
  if (!is.numeric(tau)) tau = tau_mat
  #AVE
  AVE   = list(AVE_X           = AVE_X, 
               AVE_outer_model = AVE_outer,
               AVE_inner_model = AVE_inner)
  
  #output
  out = list(Y      = Y,
              a      = a, 
              b      = b, 
              c      = c, 
              astar  = astar,
              C      = C, 
              tau    = tau, 
              scheme = scheme,
              ncomp  = ncomp, 
              crit   = crit,
              mode   = mode,
              AVE    = AVE, 
              A_mean = A_mean, 
              A_sd   = A_sd)

  class(out) = "mgcca"
  return(out)
}