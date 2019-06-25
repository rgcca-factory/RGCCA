ci.mgcca = function(object, A, BB, alpha = 0.05, ndim = 1, verbose = FALSE, plot = FALSE, nb_cores = 4, paired = F){
  
  n    = NROW(A[[1]])
  L    = length(A)
  DIM  = lapply(A, dim)
  LEN  = unlist(lapply(DIM, length))
  B_3D = which(LEN == 3)
  
  
  ##############################################
  # Initialization of the outer weight vectors #
  ##############################################
  if(any(object$ncomp!=1)){
    Yinit <- sapply(object$Y, function(x) x[, ndim])
  }else{
    Yinit <- sapply(object$Y, cbind)
  } 
  
  W = lapply(object$a, function(x) x[, ndim])
  B = lapply(object$b[B_3D], function(x) x[, ndim])
  C = lapply(object$c[B_3D], function(x) x[, ndim])
  
  ########################################
  # Construction of the Boostrap samples #
  ########################################
  if (verbose){
    boot_b  = pbmcapply::pbmclapply(1:BB, function(z) m_bootstrap(n = n, L = L, A = A, object = object, 
                                                                  W = W, B = B, C = C, 
                                                                  ndim = ndim, B_3D = B_3D, paired = paired), 
                                    mc.cores = nb_cores)
  }else{
    boot_b  = parallel::mclapply(1:BB, function(z) m_bootstrap(n = n, L = L, A = A, object = object, 
                                                               W = W, B = B, C = C, 
                                                               ndim = ndim, B_3D = B_3D, paired = paired), 
                                 mc.cores = nb_cores)
  }
  
  ########################################
  # Construction of the Coefficeints'CI  #
  ########################################

  Astar = NULL
  for (i in 1:length(boot_b[[1]]$Astar)){
    Astar[[i]] = sapply(1:BB, FUN = function(x) boot_b[[x]]$Astar[[i]])
  }
  
  if (class(Astar) == "numeric"){
    M1_Astar = c(mean(Astar), sd(Astar))
  }else{
    M1_Astar = lapply(Astar, function(w) ifelse(class(w) == "numeric", 
                                                return(matrix(c(mean(w), sd(w)), nrow = 2)), 
                                                return(apply(w, 1,  function(x) c(mean(x), sd(x))))))
  }
  
  mat_Astar  = list()
  tail       = qnorm(1-alpha/(2))
  for (j in 1:length(boot_b[[1]]$Astar)){
    mat_Astar[[j]]     <- cbind(W[[j]], M1_Astar[[j]][1, ]-tail*M1_Astar[[j]][2, ], 
                                M1_Astar[[j]][1, ]+tail*M1_Astar[[j]][2, ])
    rownames(mat_Astar[[j]]) <- colnames(A[[j]])
    colnames(mat_Astar[[j]]) <- c("Initial weights","Lower Bound","Upper Bound")
  }
  
  Bstar = NULL
  Cstar = NULL
  for (i in 1:length(boot_b[[1]]$Bstar)){
    Bstar[[i]] = sapply(1:BB, FUN = function(x) boot_b[[x]]$Bstar[[i]])
    Cstar[[i]] = sapply(1:BB, FUN = function(x) boot_b[[x]]$Cstar[[i]])
  }
  
  
  M1_Bstar = lapply(Bstar, function(w) apply(w, 1,  function(x) c(mean(x), sd(x))))
  M1_Cstar = lapply(Cstar, function(w) apply(w, 1,  function(x) c(mean(x), sd(x))))
  
  mat_Bstar  = list()
  mat_Cstar  = list()
  for (j in 1:length(boot_b[[1]]$Bstar)){
    mat_Bstar[[j]]     <- cbind(B[[j]], M1_Bstar[[j]][1, ]-tail*M1_Bstar[[j]][2, ], 
                                M1_Bstar[[j]][1, ]+tail*M1_Bstar[[j]][2, ])
    mat_Cstar[[j]]     <- cbind(C[[j]], M1_Cstar[[j]][1, ]-tail*M1_Cstar[[j]][2, ], 
                                M1_Cstar[[j]][1, ]+tail*M1_Cstar[[j]][2, ])
    rownames(mat_Bstar[[j]]) <- colnames(A[[j]])
    rownames(mat_Cstar[[j]]) <- colnames(A[[j]])
    colnames(mat_Bstar[[j]]) <- c("Initial weights","Lower Bound","Upper Bound")
    colnames(mat_Cstar[[j]]) <- c("Initial weights","Lower Bound","Upper Bound")
  }
  

  ########################################
  # Construction of the Correlations'CI  #
  ########################################
  MAT_COR = t(sapply(1:BB, FUN = function(x) boot_b[[x]][[1]]))
  M2      = apply(MAT_COR, 2,  function(x) c(mean(x), sd(x)))
  
  inner_relation           = matrix(0, 3, (L*(L-1))/2)
  connection               = cbind(rep(paste0("X", 1:L), L), rep(paste0("X", 1:L), each = L))
  colnames(inner_relation) = matrix(paste0(connection[, 1], "-", connection[, 2]), L, L)[upper.tri(diag(L))]
  rownames(inner_relation) = c("Initial Coorelation","Lower Bound","Upper Bound")
  inner_relation[1, ]      = (cor(Yinit)*object$C)[upper.tri(diag(L))]
  inner_relation[2, ]      = (matrix(M2[1, ]-tail*M2[2, ], L, L)*object$C)[upper.tri(diag(L))]
  inner_relation[3, ]      = (matrix(M2[1, ]+tail*M2[2, ], L, L)*object$C)[upper.tri(diag(L))]
  inner_relation           = t(inner_relation)
  
  return(list(a_boot = Astar, a_CI = mat_Astar,
              b_boot = Bstar, b_CI = mat_Bstar,
              c_boot = Cstar, c_CI = mat_Cstar,
              inner_relation = inner_relation))
}
