#' The function sgccak() is called by sgcca() and does not have to be used by the user. 
#' sgccak() enables the computation of SGCCA block components, outer weight vectors, etc., 
#' for each block and each dimension. 
#' @param A  A list that contains the \eqn{J} blocks of variables from which block components are constructed. 
#' It could be eiher the original matrices (\eqn{X_1, X_2, ..., X_J}) or the residual matrices (\eqn{X_{h1}, X_{h2}, ..., X_{hJ}}).
#' @param C  A design matrix that describes the relationships between blocks.
#' @param c1 A \eqn{1 * J} vector that contains the value of c1 applied to each block. The L1 bound on a[[j]] is 
#' \deqn{ \|a_{j}\|_{\ell_1} \leq c_1[j] \sqrt{p_j}.}
#' with \eqn{p_j} the number of variables of \eqn{X_j} and with c1[j] between 0 and 1 (larger L1 bound corresponds to less penalization).
#' @param scheme  Either "horst", "factorial" or "centroid" (default: centroid).
#' @param scale  If scale = TRUE, each block is standardized to zero means and unit variances (default: TRUE).
#' @param init Mode of initialization of the SGCCA algorithm. Either by Singular Value Decompostion ("svd") or random ("random") (default: "svd").
#' @param bias Logical value for biaised (\eqn{1/n}) or unbiaised (\eqn{1/(n-1)}) estimator of the var/cov.
#' @param verbose  Reports progress while computing, if verbose = TRUE (default: TRUE).
#' @param tol Stopping value for convergence.
#' @return \item{Y}{A \eqn{n * J} matrix of SGCCA block components.}
#' @return \item{a}{A list of \eqn{J} elements. Each element contains the outer weight vector of each block.}
#' @return \item{crit}{The values of the objective function at each iteration of the iterative procedure.}
#' @return \item{converg}{Speed of convergence of the alogrithm to reach the tolerance.}
#' @return \item{AVE}{Indicators of model quality based on the Average Variance Explained (AVE): AVE(for one block), AVE(outer model), AVE(inner model).}
#' @return \item{C}{A design matrix that describes the relationships between blocks (user specified).}
#' @return \item{scheme}{The scheme chosen by the user (user specified).}
#' @title Internal function for computing the SGCCA parameters (SGCCA block components, outer weight vectors etc.)
#' @export sgccak 
sgccak <-  function(A, C, c1 = rep(1, length(A)), scheme = "centroid", scale = FALSE,
                    tol = .Machine$double.eps, init="svd", bias = TRUE, verbose = TRUE){

  J <- length(A)
  pjs = sapply(A, NCOL)
  AVE_X <- rep(0, J)
  # Data standardization 
  #if (scale == TRUE) A <- lapply(A, function(x) scale2(x, bias = bias))
  #  Choose J arbitrary vectors
  if (init=="svd") {
    #SVD Initialisation of a_j or \alpha_j
    a <- lapply(A, function(x) return(svd(x,nu=0,nv=1)$v)) #
  } else if (init=="random") {
    a <- lapply(pjs,rnorm)
  } else {
    stop("init should be either random or svd.")
  }

  if (any( c1 < 1/sqrt(pjs) | c1 > 1 )) stop("L1 constraints must vary between 1/sqrt(p_j) and 1.") 
  
  const <- c1*sqrt(pjs)
  #	Apply the constraints of the general otpimization problem
  #	and compute the outer components
  iter <- 1
  converg <- crit <- numeric()
  Y <- Z <- matrix(0,NROW(A[[1]]),J)
  for (q in 1:J){
      Y[,q] <- apply(A[[q]],1,miscrossprod,a[[q]])
      a[[q]] <- soft.threshold(a[[q]], const[q])   
      a[[q]] <- as.vector(a[[q]])/norm2(a[[q]])
  }	 
  a_old <- a
  
  # Compute the value of the objective function
  
  ifelse((mode(scheme) != "function"), {g <- function(x) switch(scheme,horst=x,factorial=x**2,centroid=abs(x))
  crit_old <- sum(C*g(cov2(Y, bias = bias)))}
  ,  crit_old <- sum(C*scheme(cov2(Y, bias = bias)))
  )
  
  if (mode(scheme) == "function") 
    dg = Deriv::Deriv(scheme, env = parent.frame())
  
  repeat{

      for (q in 1:J){
      
        if (mode(scheme) == "function") {
          dgx = dg(cov2(Y[, q], Y, bias = bias))
          CbyCovq = C[q, ]*dgx
        } 
        
        else{
          if (scheme == "horst"    ) CbyCovq <- C[q, ]
          if (scheme == "factorial") CbyCovq <- C[q, ]*2*cov2(Y, Y[, q], bias = bias)
          if (scheme == "centroid" ) CbyCovq <- C[q, ]*sign(cov2(Y, Y[,q], bias = bias))
        }
      
        Z[,q] <- rowSums(mapply("*", CbyCovq,as.data.frame(Y)))
        a[[q]] <- apply( t(A[[q]]),1,miscrossprod, Z[,q])
        a[[q]] <- soft.threshold(a[[q]], const[q])
        a[[q]] <- as.vector(a[[q]])/norm2(a[[q]]) 
        Y[,q] <- apply(A[[q]], 1, miscrossprod,a[[q]])
      }
    
    #check for convergence of the SGCCA alogrithm to a solution of the stationnary equations
    
    ifelse((mode(scheme) != "function"), {g <- function(x) switch(scheme,horst=x,factorial=x**2,centroid=abs(x))
    crit[iter] <- sum(C*g(cov2(Y, bias = bias)))}
    , crit[iter] <- sum(C*scheme(cov2(Y, bias = bias)))
    )
    
    # Print out intermediate fit
    
    if (verbose & (iter %% 1)==0)
      cat(" Iter: ", formatC(iter,width=3, format="d"),
          " Fit: ",  formatC(crit[iter], digits=8, width=10, format="f"),
          " Dif: ",  formatC(crit[iter]-crit_old, digits=8, width=10, format="f"),
          "\n")
    
    stopping_criteria = c(drop(crossprod(Reduce("c", mapply("-", a, a_old))))
                          , abs(crit[iter]-crit_old))
    
    if ( any(stopping_criteria < tol) | (iter > 1000)) 
      break
    
    crit_old = crit[iter]
    a_old <- a
    iter <- iter + 1 
  }
  
  
  if (iter > 1000) warning("The SGCCA algorithm did not converge after 1000 iterations.")
  if(iter<1000 & verbose) cat("The SGCCA algorithm converged to a stationary point after", iter-1, "iterations \n") 
  if (verbose) plot(crit, xlab = "iteration", ylab = "criteria")
  
  for (q in 1:J) if(sum(a[[q]]!=0) <= 1) warning(sprintf("Deflation failed because only one variable was selected for block #",q))
  
  AVE_inner  <- sum(C*cor(Y)^2/2)/(sum(C)/2) # AVE inner model
  
  result <- list(Y = Y, a = a, crit = crit[which(crit != 0)], 
                 AVE_inner = AVE_inner, C = C, c1, scheme = scheme)
  return(result)                                                                          
}


