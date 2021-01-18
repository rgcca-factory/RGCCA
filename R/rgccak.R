#' The function rgccak() is called by rgccad() and does not have to be used by 
#' the user. The function rgccak() computes the RGCCA block components, block 
#' weight vectors, etc., for each block and each dimension. Depending on the 
#' dimensionality of each block \eqn{X_j , j = 1, ..., J}, the primal 
#' (when \eqn{n > p_j}) or the dual (when \eqn{n < p_j}) algorithm is used 
#' (see Tenenhaus et al. 2015) 
#' @inheritParams select_analysis
#' @inheritParams rgccaNa
#' @inheritParams rgccad
#' @param A  A list that contains the \eqn{J} blocks of variables from which 
#' block components are constructed. It could be either the original matrices 
#' (\eqn{X_1, X_2, ..., X_J}) or the residual matrices 
#' (\eqn{X_{h1}, X_{h2}, ..., X_{hJ}}).
#' @param na.rm If TRUE, RGCCA is run only on available data (default value) 
#' otherwise the NIPALS algorithm is used.
#' @return \item{Y}{A \eqn{n * J} matrix of block components}
#' @return \item{Z}{A \eqn{n * J} matrix of inner components}
#' @return \item{a}{A list of \eqn{J} elements. Each element of \eqn{a} is a 
#' matrix that contains the block weight vectors for each block.}
#' @return \item{converg}{Convergence of the algorithm to reach the tolerance.}
#' @return \item{AVE}{A list of numerical values giving the indicators of model 
#' quality based on the Average Variance Explained (AVE): AVE(for each block), 
#' AVE(outer model), AVE(inner model).}
#' @return \item{call}{Call of the function}
#' @return \item{crit}{A list of max(ncomp) elements. Each element 
#' (one per deflation stage) is a vector that contains the value of the RGCCA 
#' objective function across iterations.}
#' @return \item{tau}{Either a 1*J vector or a \eqn{\mathrm{max}(ncomp)\times J} 
#' matrix containing the values of the regularization parameters . The shrinkage 
#' parameter tau varies from 0 (maximizing the correlation) to 1 (maximizing the 
#' covariance). If tau = "optimal" the regularization paramaters are estimated 
#' for each block and each dimension using the Schafer and Strimmer (2005)
#' analytical formula . If tau is a \eqn{1\times J} vector, tau[j] is identical 
#' across the dimensions of block \eqn{\mathbf{X}_j}. If tau is a matrix, 
#' tau[k, j] is associated with \eqn{\mathbf{X}_{jk}} (\eqn{k}th residual matrix 
#' for block \eqn{j}). It can also be estimated by using 
#' \link{rgcca_permutation}.}
#' @references Tenenhaus M., Tenenhaus A. and Groenen PJF (2017), Regularized 
#' generalized canonical correlation analysis: A framework for sequential 
#' multiblock component methods, Psychometrika, 82, 737–777
#' @references Tenenhaus A., Philippe C., and Frouin V. (2015). Kernel 
#' Generalized Canonical Correlation Analysis. Computational Statistics and 
#' Data Analysis, 90, 114-131.
#' @references Tenenhaus A. and Tenenhaus M., (2011), Regularized Generalized 
#' Canonical Correlation Analysis, Psychometrika, Vol. 76, Nr 2, pp 257-284.
#' @references Schafer J. and Strimmer K., (2005), A shrinkage approach to 
#' large-scale covariance matrix estimation and implications for functional 
#' genomics. Statist. Appl. Genet. Mol. Biol. 4:32.
#' @title Internal function for computing the RGCCA parameters (RGCCA block 
#' components, outer weight vectors, etc.).
#' @importFrom MASS ginv
#' @importFrom stats cor rnorm
#' @importFrom graphics plot
#' @importFrom Deriv Deriv

rgccak=function (A, C, tau = "optimal", scheme = "centroid", verbose = FALSE, 
                 init = "svd", bias = TRUE, tol = 1e-08, na.rm = TRUE, 
                 scale = TRUE, scale_block = TRUE 
                 )
{  
  call = list(A = A, C = C, scheme = scheme, verbose = verbose, init = init, 
              bias = bias, tol = tol, na.rm = na.rm,
              scale = scale, scale_block = scale_block)
        
  if(mode(scheme) != "function") 
  {
    if(!scheme %in% c("horst", "factorial", "centroid"))
      {stop_rgcca("Please choose scheme as 'horst', 'factorial', 'centroid'")}
    if(scheme == "horst"){ g <- function(x) x}
    if(scheme == "factorial"){ g <- function(x)  x^2}  
    if(scheme == "centroid"){g <- function(x) abs(x)}
} 
  else g <- scheme
    
    
    J <- length(A) # number of blocks
    n <- NROW(A[[1]]) # number of individuals
    pjs <- sapply(A, NCOL) # number of variables per block
    Y <- matrix(0, n, J)
    if (!is.numeric(tau)) 
        tau = sapply(A, tau.estimate) # From Schafer and Strimmer, 2005
    
    A0 <-A
    A <- lapply(A, as.matrix)
    a <- alpha <- M <- Minv <- K <- list()
    
    # Whether primal or dual for each block
    which.primal <- which((n >= pjs) == 1) 
    which.dual <- which((n < pjs) == 1)

    # Initialisation by SVD
    if (init == "svd") { 
        for (j in which.primal) {
            a[[j]] <- initsvd(A[[j]]) 
        }
        for (j in which.dual) { 
            alpha[[j]] <- initsvd(A[[j]])
            K[[j]] <- pm( A[[j]] , t(A[[j]]), na.rm = na.rm)
        }
    }
    else if (init == "random") {
        for (j in which.primal) {
            a[[j]] <- rnorm(pjs[j]) # random initialisation
        }
        for (j in which.dual) {
            alpha[[j]] <- rnorm(n)
            K[[j]] <- pm(A[[j]] , t(A[[j]]), na.rm = na.rm)
        }
    }
    else {
        stop_rgcca("init should be either random or by SVD.")
    }
   
    N = ifelse(bias, n, n - 1)
	
    for (j in which.primal) 
    {
        ifelse(tau[j] == 1,
        {
            a[[j]] <- drop(1/sqrt(t(a[[j]]) %*% a[[j]])) * a[[j]] 
            if(a[[j]][1]<0){a[[j]]=-a[[j]]}
            Y[, j] <- pm(A[[j]] , a[[j]], na.rm = na.rm)
        }, 
        {
           
            M[[j]] <- ginv(tau[j] * diag(pjs[j]) + ((1 - tau[j]))*1/N*
                             (pm(t(A[[j]]), A[[j]], na.rm = na.rm))) 
            a[[j]] <- drop(1/sqrt(t(a[[j]])%*% M[[j]]%*%a[[j]]))*
                             (M[[j]]%*%a[[j]]) 
            if(a[[j]][1]<0){a[[j]] = -a[[j]]}
            Y[, j] <-pm(A[[j]], a[[j]], na.rm = na.rm)
        })
    }
    for (j in which.dual)
    {
        ifelse(tau[j] == 1, {
            alpha[[j]] = drop(1/sqrt(t(alpha[[j]])%*%K[[j]]%*%
                                       alpha[[j]]))*alpha[[j]]
            a[[j]] = pm(t(A[[j]]), alpha[[j]], na.rm = na.rm)
            if(a[[j]][1]<0){a[[j]] = -a[[j]]}
            Y[, j] = pm(A[[j]], a[[j]], na.rm = na.rm)
        }, {
          
            M[[j]] <- tau[j] * diag(n) + ((1 - tau[j])) *1/N* K[[j]] 
            Minv[[j]] = ginv(M[[j]])
            alpha[[j]] = drop(1/sqrt(t(alpha[[j]])%*% 
                              M[[j]]%*% K[[j]]%*% alpha[[j]])) * alpha[[j]]
            a[[j]] = pm( t(A[[j]]), alpha[[j]],na.rm=na.rm) 
            if(a[[j]][1]<0){a[[j]] = -a[[j]]}
            Y[, j] = pm(A[[j]], a[[j]], na.rm=na.rm) 
        })
    }

    crit_old <- sum(C * g(cov2(Y, bias = bias)), na.rm = na.rm)
    iter = 1
    crit = numeric()
    Z = matrix(0, NROW(A[[1]]), J)
    a_old = a
    
    dg = Deriv::Deriv(g, env = parent.frame())
   
    repeat 
    { 
      Yold <- Y 

       for (j in which.primal)
      { 
         dgx = dg(cov2(Y[, j], Y, bias = bias))
         if(tau[j] == 1)
          { 
           Z[, j] = rowSums(matrix(rep(C[j, ], n), n, J, byrow = TRUE)*
                              matrix(rep(dgx, n), n, J, byrow = TRUE)*Y, 
                            na.rm = na.rm)
		       a[[j]] = drop(1/sqrt(pm(pm(t(Z[, j]), A[[j]], na.rm = na.rm),  
		                               pm(t(A[[j]]), Z[, j], na.rm = na.rm),
		                               na.rm = na.rm)))*
		                        pm(t(A[[j]]), Z[,  j], na.rm = na.rm)  
		     if(a[[j]][1]<0){a[[j]]=-a[[j]]}

		     Y[, j] = pm(A[[j]], a[[j]], na.rm = na.rm)
			       
           }else
			      { 
             Z[, j] = rowSums(matrix(rep(C[j, ], n), n,  J, byrow = TRUE)*
                                matrix(rep(dgx, n), n,  J, byrow = TRUE)* 
                                Y, na.rm = na.rm)
             a[[j]] = drop(1/sqrt(pm(pm(t(Z[, j]), A[[j]], na.rm = na.rm), 
                                     pm(pm(M[[j]], t(A[[j]]), na.rm = na.rm), 
                                        Z[, j], na.rm = na.rm), na.rm = na.rm)))*
                      pm(M[[j]], pm(t(A[[j]]), Z[, j]))
             if(a[[j]][1]<0){a[[j]] = -a[[j]]}
             Y[, j] = pm(A[[j]], a[[j]], na.rm = na.rm)
          }
       }

      for (j in which.dual)
      {
          dgx = dg(cov2(Y[, j], Y, bias = bias))
          ifelse(tau[j] == 1, 
            {
              Z[, j] = rowSums(matrix(rep(C[j, ], n), n, J, byrow = TRUE)*
                                 matrix(rep(dgx, n), n, J, byrow = TRUE)*Y,
                               na.rm = na.rm)
              alpha[[j]] = drop(1/sqrt(t(Z[, j])%*%K[[j]]%*%Z[, j]))*Z[, j]
              a[[j]] = pm(t(A[[j]]), alpha[[j]], na.rm = na.rm)
              if(a[[j]][1]<0){a[[j]] = -a[[j]]}
              Y[, j] = pm(A[[j]], a[[j]], na.rm = na.rm)
           }, 
           {
            Z[, j] = rowSums(matrix(rep(C[j, ], n), n, J, byrow = TRUE)*
                               matrix(rep(dgx, n), n, J, byrow = TRUE)*
                               Y, na.rm = na.rm)
          	alpha[[j]] = drop(1/sqrt(t(Z[, j])%*%K[[j]]%*%Minv[[j]]%*%Z[, j]))*
          	             (Minv[[j]]%*%Z[,  j])
                   
		   a[[j]] = pm(t(A[[j]]), alpha[[j]], na.rm = na.rm)
		   if(a[[j]][1]<0){a[[j]] = -a[[j]]}
            Y[, j] = pm( A[[j]], a[[j]], na.rm = na.rm)
          })
      }

      crit[iter] <- sum(C * g(cov2(Y, bias = bias)), na.rm = na.rm)
      if (verbose & (iter%%1) == 0) 
      {
          cat(" Iter: ", formatC(iter, width = 3, format = "d"), 
              " Fit:", formatC(crit[iter], digits = 8, 
                               width = 10, format = "f"), 
              " Dif: ", formatC(crit[iter] - crit_old, digits = 8, 
                                width = 10, format = "f"), "\n")
      }
       stopping_criteria = c(drop(crossprod(Reduce("c", 
                                                   mapply("-", a, a_old)))), 
                             crit[iter] - crit_old)
     
       if (any(stopping_criteria < tol) | (iter > 1000))
        {break}  
      crit_old = crit[iter]
      a_old <- a
      iter <- iter + 1
    }
    if (iter > 1000) 
        warning("The RGCCA algorithm did not converge after 1000 iterations.")
    if (iter < 1000 & verbose) 
        cat("The RGCCA algorithm converged to a stationary point after",  
            iter - 1, "iterations \n")
    if (verbose) 
    {
        plot(crit[1:iter], xlab = "iteration", ylab = "criteria")
    }
    AVEinner <- sum(C * cor(Y)^2/2)/(sum(C)/2)
    call$tau=tau
    
    result <- list(Y = Y, a = a, crit = crit, AVE_inner = AVEinner,call=call,tau=tau)
    return(result)
}