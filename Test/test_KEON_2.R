library(Matrix)
library(utils)
library(cba)
library(RGCCA)

mckeon  = function(A)
{
  J   = length(A)
  SD  = bdiag(lapply(A, cov2))
  S   = cov2(Reduce("cbind", A))
  rho = (1/(J-1))*(eigen(solve(SD)%*%S)$values[1]-1)
  return(Re(rho))
}

set.seed(123456789)
d1 = runif(50, 0, 10)
d2 = runif(50, 0, 10)
E = matrix(rnorm(800, 0, 1), 50, 16)

X1 = matrix(0, 50, 4)
X1[, 1] = d1 + E[, 1]
X1[, 2] = d1 + E[, 2]
X1[, 3] = d1 + E[, 3]
X1[, 4] = d1 + E[, 4]
colnames(X1) = paste("X1", 1:4, sep = "")

X2 = matrix(0, 50, 2)
X2[, 1] = 0.7*d1 + E[, 5]
X2[, 2] = 0.7*d1 + E[, 6]
colnames(X2) = paste("X1", 1:2, sep = "")

X3 = matrix(0, 50, 2)
X3[, 1] = 0.4*d1 + E[, 7]
X3[, 2] = 0.4*d1 + E[, 8]
colnames(X3) = paste("X3", 1:2, sep = "")

X4 = matrix(0, 50, 3)
X4[, 1] = d2 + E[, 9]
X4[, 2] = d2 + E[, 10]
X4[, 3] = d2 + E[, 11]
colnames(X4) = paste("X4", 1:3, sep = "")

X5 = matrix(0, 50, 3)
X5[, 1] = 0.7*d2 + E[, 12]
X5[, 2] = 0.7*d2 + E[, 13]
X5[, 3] = 0.7*d2 + E[, 14]
colnames(X5) = paste("X5", 1:3, sep = "")

X6 = matrix(0, 50, 2)
X6[, 1] = 0.4*d2 + E[, 15]
X6[, 2] = 0.4*d2 + E[, 16]
colnames(X6) = paste("X6", 1:2, sep = "")

A = list(X1 = X1, X2 = X2, X3 = X3,
         X4 = X4, X5 = X5, X6 = X6)



dendo_KEON_2 = function(A){
  # browser()
  J = length(A)
  
  TREE = NULL
  TREE = lapply(J:2, function(x) append(TREE, list(idx_2D                 = NULL, 
                                                   idx_2D_max_KEON        = 0, 
                                                   A_current              = NULL,
                                                   idx_melt               = 0,
                                                   idx_2D_max_KEON_first  = 0,
                                                   max_KEON               = 0,
                                                   max_KEON_first         = 0,
                                                   fixed                  = F)))
  
  out        = list()
  out$labels = names(A)

  names(A) = -(1:J)
  
  TREE[[1]]$idx_2D    = combn(x = 1:J, m = 2)
  TREE[[1]]$A_current = A
  
  i = 1

  while (i < J ){
    
    J_current = J - i + 1
    
    if (i != 1 && TREE[[i]]$idx_2D_max_KEON_first == 0){
      A_current                     = TREE[[i - 1]]$A_current
      A_current[[J_current+2]]      = cbind(A_current[[TREE[[i - 1]]$idx_melt[1]]], 
                                            A_current[[TREE[[i - 1]]$idx_melt[2]]])
      names(A_current)[J_current+2] = i - 1
      A_current[TREE[[i - 1]]$idx_melt] = NULL
    }else{
      A_current = TREE[[i]]$A_current
    }
     
    # browser()
    A_current_scaled = lapply(A_current, scale2)
    # A_current_scaled = lapply(A_current_scaled, function(x) x/sqrt(NCOL(x)))
    
    
    if (TREE[[i]]$idx_2D_max_KEON_first == 0){
      idx_2D = combn(x = 1:J_current, m = 2)
    }else{
      idx_2D = TREE[[i]]$idx_2D
    }
    KEON            = apply(X = idx_2D, 2, FUN = function(x) mckeon(A_current_scaled[x]))
    idx_2D_max_KEON = which.max(KEON)
    idx_melt        = idx_2D[ , idx_2D_max_KEON]
    max_KEON        = max(KEON)
    
    melt_1 = as.numeric(names(A_current)[idx_melt[1]])
    if (melt_1 < 0){
      max_1   = 1
      fixed_1 = F
    }else{
      max_1   = TREE[[melt_1]]$max_KEON
      fixed_1 = TREE[[melt_1]]$fixed
    }
    
    melt_2 = as.numeric(names(A_current)[idx_melt[2]])
    if (melt_2 < 0){
      max_2   = 1
      fixed_2 = F
    }else{
      max_2   = TREE[[melt_2]]$max_KEON
      fixed_2 = TREE[[melt_2]]$fixed
    }
    
    if (i != 1){
      if ( any(max_KEON > c(max_1, max_2)) || any(c(fixed_1, fixed_2)) ){
        if (dim(TREE[[i - 1]]$idx_2D)[2] != 1){
          # browser()
          TREE[[i - 1]]$idx_2D = as.matrix(TREE[[i - 1]]$idx_2D[, -TREE[[i - 1]]$idx_2D_max_KEON])

          i = i - 1
        }else{
          if (i !=2){
            DIM = unlist(lapply(1:(i-2), function(x) dim(TREE[[x]]$idx_2D)[2]))
            if (max(unique(DIM)) == 1){
              browser()
              j                         = min(unlist(lapply(1:(i-2), function(x) which(TREE[[x]]$fixed == F))))
              TREE[[j]]$idx_2D_max_KEON = TREE[[j]]$idx_2D_max_KEON_first
              TREE[[j]]$max_KEON        = TREE[[j]]$max_KEON_first
              idx_2D                    = combn(x = 1:(J_current + j), m = 2)
              TREE[[j]]$idx_melt        = idx_2D[, TREE[[j]]$idx_2D_max_KEON_first]
              TREE[[j]]$fixed           = T
              
              TREE[(j + 1) : (J - 1)] = lapply((j + 1) : (J - 1), function(x) list(idx_2D                 = NULL, 
                                                                                    idx_2D_max_KEON        = 0, 
                                                                                    A_current              = NULL,
                                                                                    idx_melt               = 0,
                                                                                    idx_2D_max_KEON_first  = 0,
                                                                                    max_KEON               = 0,
                                                                                    max_KEON_first         = 0,
                                                                                    fixed                  = F))
              
              i = j + 1
            }else{
              browser()
              j                = max(which(DIM != 1))
              TREE[[j]]$idx_2D = as.matrix(TREE[[j]]$idx_2D[, -TREE[[j]]$idx_2D_max_KEON])
              
              TREE[(j + 1) : (J - 1)] = lapply((j + 1) : (J - 1), function(x) list(idx_2D                 = NULL, 
                                                                                    idx_2D_max_KEON        = 0, 
                                                                                    A_current              = NULL,
                                                                                    idx_melt               = 0,
                                                                                    idx_2D_max_KEON_first  = 0,
                                                                                    max_KEON               = 0,
                                                                                    max_KEON_first         = 0,
                                                                                    fixed                  = F))
              
              i = j
            }
          }else{
            if (dim(TREE[[1]]$idx_2D)[2] == 1){
              TREE[[1]]$idx_2D_max_KEON = TREE[[1]]$idx_2D_max_KEON_first
              TREE[[1]]$max_KEON        = TREE[[1]]$max_KEON_first
              idx_2D                    = combn(x = 1:J, m = 2)
              TREE[[1]]$idx_melt        = idx_2D[, TREE[[1]]$idx_2D_max_KEON_first]
              TREE[[1]]$fixed           = T
              
              TREE[[2]] = list(idx_2D                 = NULL, 
                             idx_2D_max_KEON        = 0, 
                             A_current              = NULL,
                             idx_melt               = 0,
                             idx_2D_max_KEON_first  = 0,
                             max_KEON               = 0,
                             max_KEON_first         = 0,
                             fixed                  = F)
              
              i = 2
            }else{
              TREE[[1]]$idx_2D = as.matrix(TREE[[1]]$idx_2D[, -TREE[[1]]$idx_2D_max_KEON])
              
              TREE[[2]] = list(idx_2D                 = NULL, 
                             idx_2D_max_KEON        = 0, 
                             A_current              = NULL,
                             idx_melt               = 0,
                             idx_2D_max_KEON_first  = 0,
                             max_KEON               = 0,
                             max_KEON_first         = 0,
                             fixed                  = F)
              
              i = 1
            }
          }
        }
      }else{
        TREE[[i]]$idx_2D_max_KEON = idx_2D_max_KEON
        TREE[[i]]$idx_melt        = idx_melt
        TREE[[i]]$max_KEON        = max_KEON
        
        if (TREE[[i]]$idx_2D_max_KEON_first == 0){
          TREE[[i]]$idx_2D_max_KEON_first  = idx_2D_max_KEON
          TREE[[i]]$max_KEON_first         = max_KEON
          TREE[[i]]$A_current              = A_current
          TREE[[i]]$idx_2D                 = idx_2D
        }

        i = i + 1
      }
    }else{
      TREE[[1]]$idx_2D_max_KEON = idx_2D_max_KEON
      TREE[[1]]$idx_melt        = idx_melt
      TREE[[1]]$max_KEON        = max_KEON
      
      if (TREE[[1]]$idx_2D_max_KEON_first == 0){
        TREE[[1]]$idx_2D_max_KEON_first  = idx_2D_max_KEON
        TREE[[1]]$max_KEON_first         = max_KEON
        TREE[[1]]$idx_2D                 = idx_2D
      }
      
      i = i + 1
    }
    print(i)
  }
  
  out$method         = "McKeon"
  attr(out, "class") = "hclust"
  # browser()
  out$merge          = lapply(TREE, function(x) c(as.numeric(names(x$A_current)[x$idx_melt[1]]), 
                                          as.numeric(names(x$A_current)[x$idx_melt[2]])))
  out$merge          = Reduce("rbind", out$merge)
  out$height         = lapply(TREE, function(x) 1 - x$max_KEON)
  
  A                                           = lapply(A, scale2)
  A                                           = lapply(A, function(x) x/sqrt(NCOL(x)))
  KEON_1                                      = combn(x = 1:J, m = 2, FUN = function(x) mckeon(A[x]))
  KEON_1_TRI                                  = matrix(0, J, J)
  KEON_1_TRI[upper.tri(KEON_1_TRI, diag = F)] = KEON_1
  KEON_1_SYM                                  = diag(J) + KEON_1_TRI + t(KEON_1_TRI)
  
  out$order  = rep(0, J)
  out$order  = order.optimal(dist = as.dist(KEON_1_SYM), merge = out$merge)$order
  
  fixed = any(lapply(TREE, function(x) x$fixed))
  
  return(list(out = out, fixed = fixed))
}

res = dendo_KEON_2(A)

plot(res$out)
