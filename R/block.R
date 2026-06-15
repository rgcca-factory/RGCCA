### Create classes
new_block <- function(x, j, na.rm = TRUE, bias = TRUE,
                      ..., class = character()) {
  n <- NROW(x)
  p <- NCOL(x)
  N <- ifelse(bias, n, n - 1)

  x <- list(
    x = x,
    j = j,
    n = n,
    p = p,
    N = N,
    na.rm = na.rm,
    a = NULL,
    Y = NULL,
    ...
  )
  class(x) <- c(class, "block")
  
  return(x)
}

new_dual_block <- function(x, j, na.rm = TRUE, ..., class = character()) {
  K <- pm(x, t(x), na.rm = na.rm)
  new_block(
    x, j, na.rm, alpha = NULL, K = K, ..., class = c(class, "dual_block")
  )
}

new_primal_regularized_block <- function(x, j, tau, ...) {
  new_block(x, j, tau = tau, M = NULL, ..., class = "primal_regularized_block")
}

new_dual_regularized_block <- function(x, j, tau, ...) {
  new_dual_block(
    x, j, tau = tau, M = NULL, ..., class = "dual_regularized_block"
  )
}

new_sparse_block <- function(x, j, sparsity, tol = 1e-08, ...) {
  
  const <- sqrt(NCOL(x)) * sparsity[[1]]
  new_block(
    x, j, sparsity = sparsity, const = const,
    tol = tol, ..., class = "sparse_block"
  )
}

new_tensor_block <- function(x, j, rank, mode_orth, ..., class = character()) {
  
  new_block(
    x, j, rank = rank, mode_orth = mode_orth, factors = NULL,
    lambda = NULL, ..., class = c(class, "tensor_block")
  )
}
new_sparse_tensor_block <- function(x, j, rank, mode_orth,sparsity,sparse_lambda,tol = 1e-08, ..., class = character()) {
 
  const<-c(rep(NA,length(dim(x))-1) )
  
 
  for (m in 1:(length(dim(x))-1)){
   
   
          const[m]<-sqrt(dim(x)[1+m]) * sparsity[m]

    
   
  }

  const2 <- sqrt(rank *sparse_lambda )
  
  
 
  new_block(
    x, j, rank = rank, sparsity=sparsity,const = const,
    tol = tol,mode_orth = mode_orth, factors = NULL,const2=const2,
    lambda = NULL, ..., class ="sparse_tensor_block"
  )
}

new_regularized_tensor_block <- function(x, j, rank, mode_orth, tau, ...) {
  new_tensor_block(
    x, j, rank = rank, mode_orth = mode_orth, tau = tau,
    M = NULL, ..., class = "regularized_tensor_block"
  )
}

new_separable_regularized_tensor_block <- function(x, j, rank, mode_orth,
                                                   tau, ...) {
  new_tensor_block(
    x, j, rank = rank, mode_orth = mode_orth, tau = tau,
    M = NULL, ..., class = "separable_regularized_tensor_block"
  )
}

### Utility method to choose the adequate class
create_block <- function(x, j, bias, na.rm, tau, sparsity,sparse_lambda,
                         tol, rank, mode_orth, separable) {
  
 
  if (length(dim(x)) > 2) {        # TGCCA
    if (tau < 1) {
      if (separable) {
        res <- new_separable_regularized_tensor_block(
          x, j, rank, mode_orth, tau, bias = bias, na.rm = na.rm
        )
      } else {
        res <- new_regularized_tensor_block(
          x, j, rank, mode_orth, tau, bias = bias, na.rm = na.rm
        )
      }
      }else {

        if  (any(unlist(sparsity)!=1)) {  
        

        
        res <- new_sparse_tensor_block(x, j,rank,  mode_orth,unlist(sparsity),sparse_lambda,tol,bias = bias, na.rm = na.rm)
        } else {
          res <- new_tensor_block(x, j, rank, mode_orth, bias = bias, na.rm = na.rm)
    }  }
    }else {
    if (sparsity < 1 &sparsity>0) {             # SGCCA
      res <- new_sparse_block(x, j, sparsity, tol, bias = bias, na.rm = na.rm)
    } else if (NROW(x) > NCOL(x)) { # Primal RGCCA
      if (tau < 1) {
        res <-
          new_primal_regularized_block(x, j, tau, bias = bias, na.rm = na.rm)
      } else {
        res <- new_block(x, j, bias = bias, na.rm = na.rm)
      }
    } else {                        # Dual RGCCA
      if (tau < 1) {
        res <- new_dual_regularized_block(x, j, tau, bias = bias, na.rm = na.rm)
      } else {
        res <- new_dual_block(x, j, bias = bias, na.rm = na.rm)
      }
    }
  }
  return(res)
}
