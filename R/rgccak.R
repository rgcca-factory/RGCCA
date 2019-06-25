#' The function rgccak() is called by rgcca() and does not have to be used by the user. 
#' The function rgccak() computes the RGCCA block components, outer weight vectors, etc., 
#' for each block and each dimension. Depending on the dimensionality of each block \eqn{X_j , j = 1, ..., J}, 
#' the primal (when \eqn{n > p_j}) or the dual (when \eqn{n < p_j}) algorithm is used (see Tenenhaus et al. 2015) 
#' @param A  A list that contains the \eqn{J} blocks of variables. Either the blocks (\eqn{X_1, X_2, ..., X_J}) or the residual matrices (\eqn{X_{h1}, X_{h2}, ..., X_{hJ}}).
#' @param C  A design matrix that describes the relationships between blocks. (Default: complete design).
#' @param tau A \eqn{1 * J} vector that contains the values of the shrinkage parameters \eqn{\tau_j}, \eqn{ j=1, ..., J}. (Default: \eqn{\tau_j = 1}, \eqn{ j=1, ..., J}).
#' If tau = "optimal" the shrinkage intensity paramaters are estimated using the Schafer and Strimmer (2005) 
#' analytical formula. 
#' @param scheme The value is "horst", "factorial", "centroid" or any diffentiable convex scheme function g designed by the user (default: "centroid").
#' @param scale  if scale = TRUE, each block is standardized to zero means and unit variances (default: TRUE).
#' @param verbose  Will report progress while computing if verbose = TRUE (default: TRUE).
#' @param init The mode of initialization to use in the RGCCA algorithm. The alternatives are either by Singular Value Decompostion or random (default : "svd").
#' @param bias A logical value for either a biaised or unbiaised estimator of the var/cov.
#' @param tol Stopping value for convergence.
#' @return \item{Y}{A \eqn{n * J} matrix of RGCCA outer components}
#' @return \item{Z}{A \eqn{n * J} matrix of RGCCA inner components}
#' @return \item{a}{A list of outer weight vectors}
#' @return \item{crit}{The values of the objective function to be optimized in each iteration of the iterative procedure.}
# #' @return \item{converg}{Speed of convergence of the algorithm to reach the tolerance.}
#' @return \item{AVE}{Indicators of model quality based on the Average Variance Explained (AVE): 
#' AVE(for one block), AVE(outer model), AVE(inner model).}
#' @return \item{C}{A design matrix that describes the relationships between blocks (user specified).}
#' @return \item{tau}{\eqn{1 * J} vector containing the value for the tau penalties applied to each of the \eqn{J} blocks of data (user specified)}
#' @return \item{scheme}{The scheme chosen by the user (user specified).}
#' @references Tenenhaus M., Tenenhaus A. and Groenen PJF (2017), Regularized generalized canonical correlation analysis: A framework for sequential multiblock component methods, Psychometrika, in press
#' @references Tenenhaus A., Philippe C., & Frouin V. (2015). Kernel Generalized Canonical Correlation Analysis. Computational Statistics and Data Analysis, 90, 114-131.
#' @references Tenenhaus A. and Tenenhaus M., (2011), Regularized Generalized Canonical Correlation Analysis, Psychometrika, Vol. 76, Nr 2, pp 257-284.
#' @references Schafer J. and Strimmer K., (2005), A shrinkage approach to large-scale covariance matrix estimation and implications for functional genomics. Statist. Appl. Genet. Mol. Biol. 4:32.
#' @title Internal function for computing the RGCCA parameters (RGCCA block components, outer weight vectors, etc.).
#' @export rgccak
#' @importFrom MASS ginv
#' @importFrom stats cor rnorm
#' @importFrom graphics plot

rgccak <- function (A, C, tau = "optimal", scheme = "centroid", scale = FALSE, 
                    verbose = FALSE, init="svd", bias = TRUE, tol = 1e-8) {
  A <- lapply(A, as.matrix)
  J <- length(A)
  n <- NROW(A[[1]])
  pjs <- sapply(A,NCOL)
  Y <- matrix(0, n, J)
  
  #if (scale == TRUE) A <- lapply(A, function(x) scale2(x, bias = bias))
  if (!is.numeric(tau)) tau = sapply(A, tau.estimate)

  
  a <- alpha <- M <- Minv <- K <- list()
  which.primal <- which((n>=pjs)==1)
  which.dual <- which((n<pjs)==1)
    
    
  ####################################
  # Initialisation and normalisation #
  ####################################
  
  if (init=="svd") {
    #SVD Initialisation of a_j or \alpha_j
    for (j in which.primal){
      a[[j]] <- initsvd(A[[j]])
    }
    for (j in which.dual){
      alpha[[j]] <- initsvd(A[[j]])
      K[[j]] <- A[[j]]%*%t(A[[j]])
    }
  } else if (init=="random") {
    #Random Initialisation of a_j or \alpha_j
    for (j in which.primal) {
      a[[j]] <- rnorm(pjs[j])
    }
    for (j in which.dual) {
      alpha[[j]] <- rnorm(n)
      K[[j]] <- A[[j]]%*%t(A[[j]])
    }
  } else {
    stop("init should be either random or by SVD.")
  }
  
  N = ifelse(bias, n, n-1)  
    
  # Normalisation of a_j or \alpha_j
  for (j in which.primal) {
    ifelse( tau[j] == 1, {
      a[[j]] <- drop(1/sqrt(t(a[[j]]) %*% a[[j]])) * a[[j]]
      Y[, j] <- A[[j]] %*% a[[j]] } , {
      M[[j]] <- ginv(tau[j] * diag(pjs[j]) + ((1 - tau[j])/(N)) * (t(A[[j]]) %*% A[[j]]))
      a[[j]] <- drop(1/sqrt(t(a[[j]]) %*% M[[j]] %*% a[[j]])) * M[[j]] %*% a[[j]]
      Y[, j] <- A[[j]] %*% a[[j]]}
    )
  }
    
  for (j in which.dual) {
      ifelse(tau[j] == 1, {alpha[[j]] = drop(1/sqrt(t(alpha[[j]])%*%K[[j]]%*% alpha[[j]]))*alpha[[j]]
                           a[[j]] = t(A[[j]])%*%alpha[[j]]
                           Y[, j] = A[[j]] %*% a[[j]]}
                        , {M[[j]] = tau[j]*diag(n)+(1-tau[j])/(N)*K[[j]]
                          Minv[[j]] = ginv(M[[j]])
                          alpha[[j]] = drop(1/sqrt(t(alpha[[j]]) %*% M[[j]]%*%K[[j]]%*% alpha[[j]]))*alpha[[j]]
                          a[[j]] = t(A[[j]])%*%alpha[[j]]
                          Y[, j] = A[[j]] %*% a[[j]]}
            )
  }
  
  ifelse((mode(scheme) != "function"), {h <- function(x) switch(scheme,horst=x,factorial=x**2,centroid=abs(x))
                                        crit_old <- sum(C*h(cov2(Y, bias = bias)))}
                                     ,  crit_old <- sum(C*scheme(cov2(Y, bias = bias)))
  )
    
  iter = 1
  crit = numeric()
  Z = matrix(0, NROW(A[[1]]), J)
  a_old = a
  if (mode(scheme) == "function") 
    dg = Deriv::Deriv(scheme, env = parent.frame())
  
  repeat {
    Yold <- Y
    
    if (mode(scheme) == "function"){
      
      ############
      # g scheme #
      ############
      
        for (j in which.primal){
          dgx = dg(cov2(Y[, j], Y, bias = bias))
          # assign(formalArgs(scheme), cov2(Y[, j], Y, bias = bias))
          # dgx = as.vector(attr(eval(dg), "grad"))
          ifelse(tau[j]==1, {Z[, j] = rowSums(matrix(rep(C[j, ], n), n, J, byrow = TRUE) * matrix(rep(dgx, n), n, J, byrow = TRUE) * Y)
                             a[[j]] = drop(1/sqrt(t(Z[, j]) %*% A[[j]] %*% t(A[[j]]) %*% Z[, j])) * (t(A[[j]]) %*% Z[, j])
                             Y[, j] = A[[j]] %*% a[[j]]}
                 , {Z[, j] = rowSums(matrix(rep(C[j, ], n), n, J, byrow = TRUE) * matrix(rep(dgx, n), n, J, byrow = TRUE) * Y)
                    a[[j]] = drop(1/sqrt(t(Z[, j]) %*% A[[j]] %*% M[[j]] %*% t(A[[j]]) %*% Z[, j])) * (M[[j]] %*% t(A[[j]]) %*% Z[, j])
                    Y[, j] = A[[j]] %*% a[[j]]}
          )
        }
      
        for (j in which.dual) {
          dgx = dg(cov2(Y[, j], Y, bias = bias))
          # assign(formalArgs(scheme), cov2(Y[, j], Y, bias = bias))
          # dgx = as.vector(attr(eval(dg), "grad"))
          ifelse(tau[j]==1, {Z[, j] = rowSums(matrix(rep(C[j, ], n), n, J, byrow = TRUE) * matrix(rep(dgx, n), n, J, byrow = TRUE) * Y)
                             alpha[[j]] = drop(1/sqrt(t(Z[, j]) %*% K[[j]] %*% Z[, j])) * Z[, j]
                             a[[j]] = t(A[[j]])%*% alpha[[j]]
                             Y[, j] = A[[j]] %*% a[[j]]}
                 , {Z[, j] = rowSums(matrix(rep(C[j, ], n), n, J, byrow = TRUE) * matrix(rep(dgx, n), n, J, byrow = TRUE) * Y)
                    alpha[[j]] = drop(1/sqrt(t(Z[, j]) %*% K[[j]] %*% Minv[[j]] %*% Z[, j])) * (Minv[[j]] %*% Z[, j])
                    a[[j]] = t(A[[j]])%*% alpha[[j]]
                    Y[, j] = A[[j]] %*% a[[j]]}
          )
        }
      }  
        
    else{
    ################
    # Horst Scheme #
    ################
    if (scheme == "horst") {
                                
        for (j in which.primal) {
          ifelse(tau[j]==1, {Z[, j] = rowSums(matrix(rep(C[j, ], n), n, J, byrow = TRUE) * Y)
                             a[[j]] = drop(1/sqrt(t(Z[, j]) %*% A[[j]] %*% t(A[[j]]) %*% Z[, j])) * t(A[[j]]) %*% Z[, j]
                             Y[, j] = A[[j]] %*% a[[j]]}
                          , {Z[, j] = rowSums(matrix(rep(C[j, ], n), n, J, byrow = TRUE) * Y)
                             a[[j]] = drop(1/sqrt(t(Z[, j]) %*% A[[j]] %*% M[[j]] %*% t(A[[j]]) %*% Z[, j])) * M[[j]] %*% t(A[[j]]) %*% Z[, j]
                             Y[, j] = A[[j]] %*% a[[j]]}
                )
        }
             
              
        for (j in which.dual) {
          ifelse(tau[j]==1, {Z[, j] = rowSums(matrix(rep(C[j, ], n), n, J, byrow = TRUE) * Y) 
                             alpha[[j]] = drop(1/sqrt(t(Z[, j]) %*% K[[j]] %*% Z[, j]))* Z[, j]
                             a[[j]] = t(A[[j]])%*% alpha[[j]]
                             Y[, j] = A[[j]] %*% a[[j]]}
          
                          , {Z[, j] = rowSums(matrix(rep(C[j, ], n), n, J, byrow = TRUE) * Y) 
                             alpha[[j]] = drop(1/sqrt(t(Z[, j]) %*% K[[j]] %*% Minv[[j]] %*% Z[, j]))*(Minv[[j]] %*% Z[, j])
                             a[[j]] = t(A[[j]])%*% alpha[[j]]
                             Y[, j] = A[[j]] %*% a[[j]]}
                )
        }
    }

    ####################
    # Factorial Scheme #
    ####################
    
    if (scheme == "factorial") {
                                              
        for (j in which.primal){
          ifelse(tau[j]==1, {Z[, j] = rowSums(matrix(rep(C[j, ], n), n, J, byrow = TRUE) * matrix(rep(cov2(Y[, j], Y, bias = bias), n), n, J, byrow = TRUE) * Y)
                             a[[j]] = drop(1/sqrt(t(Z[, j]) %*% A[[j]] %*% t(A[[j]]) %*% Z[, j])) * (t(A[[j]]) %*% Z[, j])
                             Y[, j] = A[[j]] %*% a[[j]]}
                          , {Z[, j] = rowSums(matrix(rep(C[j, ], n), n, J, byrow = TRUE) * matrix(rep(cov2(Y[, j], Y, bias = bias), n), n, J, byrow = TRUE) * Y)
                             a[[j]] = drop(1/sqrt(t(Z[, j]) %*% A[[j]] %*% M[[j]] %*% t(A[[j]]) %*% Z[, j])) * (M[[j]] %*% t(A[[j]]) %*% Z[, j])
                             Y[, j] = A[[j]] %*% a[[j]]}
                )
        }
          
        for (j in which.dual) {
          ifelse(tau[j]==1, {Z[, j] = rowSums(matrix(rep(C[j, ], n), n, J, byrow = TRUE) * matrix(rep(cov2(Y[, j], Y, bias = bias), n), n, J, byrow = TRUE) * Y)
                             alpha[[j]] = drop(1/sqrt(t(Z[, j]) %*% K[[j]] %*% Z[, j])) * Z[, j]
                             a[[j]] = t(A[[j]])%*% alpha[[j]]
                             Y[, j] = A[[j]] %*% a[[j]]}
                          , {Z[, j] = rowSums(matrix(rep(C[j, ], n), n, J, byrow = TRUE) * matrix(rep(cov2(Y[, j], Y, bias = bias), n), n, J, byrow = TRUE) * Y)
                            alpha[[j]] = drop(1/sqrt(t(Z[, j]) %*% K[[j]] %*% Minv[[j]] %*% Z[, j])) * (Minv[[j]] %*% Z[, j])
                            a[[j]] = t(A[[j]])%*% alpha[[j]]
                            Y[, j] = A[[j]] %*% a[[j]]}
                )
        }
    }

    ###################
    # Centroid Scheme #
    ###################
    
    if (scheme == "centroid") {
      for (j in which.primal) {
        ifelse(tau[j]==1, {Z[, j] = rowSums(matrix(rep(C[j, ], n), n, J, byrow = TRUE) * sign(matrix(rep(cov2(Y[, j], Y, bias = bias), n), n, J, byrow = TRUE)) * Y)
                           a[[j]] = drop(1/sqrt(t(Z[, j]) %*% A[[j]] %*% t(A[[j]]) %*% Z[, j])) * (t(A[[j]]) %*% Z[, j])
                           Y[, j] = A[[j]] %*% a[[j]]}
                        , {Z[, j] = rowSums(matrix(rep(C[j, ], n), n, J, byrow = TRUE) * sign(matrix(rep(cov2(Y[, j], Y, bias = bias), n), n, J, byrow = TRUE)) * Y)
                           a[[j]] = drop(1/sqrt(t(Z[, j]) %*% A[[j]] %*% M[[j]] %*% t(A[[j]]) %*% Z[, j])) * (M[[j]] %*% t(A[[j]]) %*% Z[, j])
                           Y[, j] = A[[j]] %*% a[[j]]}
              )
      }
      
      for (j in which.dual) {
        ifelse(tau[j]==1, {Z[, j] = rowSums(matrix(rep(C[j, ], n), n, J, byrow = TRUE) * sign(matrix(rep(cov2(Y[, j], Y, bias = bias), n), n, J, byrow = TRUE)) * Y)
                           alpha[[j]] = drop(1/sqrt(t(Z[, j]) %*% K[[j]] %*% Z[, j])) * Z[, j]
                           a[[j]] = t(A[[j]])%*% alpha[[j]]
                           Y[, j] = A[[j]] %*% a[[j]]}
                        , {Z[, j] = rowSums(matrix(rep(C[j, ], n), n, J, byrow = TRUE) * sign(matrix(rep(cov2(Y[, j], Y, bias = bias), n), n, J, byrow = TRUE)) * Y)
                           alpha[[j]] = drop(1/sqrt(t(Z[, j]) %*% K[[j]] %*% Minv[[j]] %*% Z[, j])) * (Minv[[j]] %*% Z[, j])
                           a[[j]] = t(A[[j]])%*% alpha[[j]]
                           Y[, j] = A[[j]] %*% a[[j]]}
        )
      }
    }
    }

    
    ifelse((mode(scheme) != "function"), {g <- function(x) switch(scheme,horst=x,factorial=x**2,centroid=abs(x))
    crit[iter] <- sum(C*g(cov2(Y, bias = bias)))}
    , crit[iter] <- sum(C*scheme(cov2(Y, bias = bias)))
    )    
    
    if (verbose & (iter %% 1)==0)
      cat(" Iter: ",formatC(iter,width=3, format="d"),
          " Fit:",  formatC(crit[iter], digits=8, width=10, format="f"),
          " Dif: ", formatC(crit[iter]-crit_old, digits=8, width=10, format="f"),
          "\n")
    
    stopping_criteria = c(drop(crossprod(Reduce("c", mapply("-", a, a_old))))
                          , crit[iter]-crit_old)
    
    if ( any(stopping_criteria < tol) | (iter > 1000)) 
      break
    
    crit_old = crit[iter]
    a_old <- a
    iter <- iter + 1
    
  }
  
  if (iter > 1000) warning("The RGCCA algorithm did not converge after 1000 iterations.")
  if(iter<1000 & verbose) cat("The RGCCA algorithm converged to a stationary point after", iter-1, "iterations \n")
  if (verbose) plot(crit[1:iter], xlab = "iteration", ylab = "criteria")
  
  AVEinner <- sum(C * cor(Y)^2/2)/(sum(C)/2)
  result <- list(Y = Y, a = a, crit = crit, 
                 AVE_inner = AVEinner, C = C, tau = tau, scheme = scheme)
  
  return(result)
}
