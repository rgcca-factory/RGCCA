#' The function sgccak() is called by sgcca() and does not have to be used by
#' the user. sgccak() enables the computation of SGCCA block components, outer
#' weight vectors, etc., for each block and each deflation stage.
#' @inheritParams select_analysis
#' @inheritParams rgccad
#' @inheritParams sgcca
#' @inheritParams rgccak
#' @return \item{Y}{A list of \eqn{J} elements. Each element of \eqn{Y} is a
#' matrix that contains the analysis components for the corresponding block.}
#' @return \item{a}{A list of \eqn{J} elements. Each element of \eqn{a} is a
#' matrix that contains the outer weight vectors for each block.}
#' @return \item{crit}{A vector of integer that contains for each component
#' the values of the analysis criteria across iterations.}
#' @return \item{converg}{Speed of convergence of the alogrithm to reach the
#' tolerance.}
#' @return \item{AVE}{A list of numerical values giving the indicators of model
#' quality based on the Average Variance Explained (AVE): AVE(for each block),
#' AVE(outer model), AVE(inner model).}
#' @title Internal function for computing the SGCCA parameters (SGCCA block
#' components, outer weight vectors etc.)
#' @importFrom Deriv Deriv
sgccak <-  function(A, C, sparsity = rep(1, length(A)), scheme = "centroid",
                    tol = .Machine$double.eps,
                    init = "svd", bias = TRUE, verbose = TRUE,
                    quiet = FALSE, na.rm = TRUE){

  J <- length(A)
  pjs = sapply(A, NCOL)

  #  Choose J arbitrary vectors
  if (init=="svd") {
    #SVD Initialisation for a_j
    a <- lapply(A, function(x) return(initsvd(x, dual = FALSE)))
    a <- lapply(a,function(x) return(as.vector(x)))
   } else if (init=="random")
    {
    a <- lapply(pjs,rnorm)
    } else {
    stop_rgcca("init should be either random or svd.")
    }

  if (any( sparsity < 1/sqrt(pjs) | sparsity > 1 ))
    stop_rgcca("L1 constraints must vary between 1/sqrt(p_j) and 1.")

  const <- sparsity*sqrt(pjs)
  #	Apply the constraints of the general otpimization problem
  #	and compute the outer components
  iter <- 1
   crit <- numeric()
  Y <- Z <- matrix(0,NROW(A[[1]]),J)
  for (q in 1:J){
      a[[q]] <- soft.threshold(a[[q]], const[q])
      a[[q]] <- as.vector(a[[q]])/norm2(a[[q]])
      Y[, q] <- pm(A[[q]], a[[q]], na.rm = na.rm)
  }
  a_old <- a

  # Compute the value of the objective function
  ifelse((mode(scheme) != "function"),
         {g <- function(x) switch(scheme,
                                  horst = x,
                                  factorial = x**2,
                                  centroid = abs(x))
          crit_old <- sum(C*g(cov2(Y, bias = bias)))
          },
         crit_old <- sum(C*scheme(cov2(Y, bias = bias)))
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
          if (scheme == "horst")
            CbyCovq <- C[q, ]
          if (scheme == "factorial")
            CbyCovq <- C[q, ]*2*cov2(Y, Y[, q], bias = bias)
          if (scheme == "centroid")
            CbyCovq <- C[q, ]*sign(cov2(Y, Y[,q], bias = bias))
        }

        Z[, q] <- rowSums(mapply("*", CbyCovq,as.data.frame(Y)))
        a[[q]] <- pm(t(A[[q]]), Z[, q], na.rm = na.rm)
        a[[q]] <- soft.threshold(a[[q]], const[q])
        a[[q]] <- as.vector(a[[q]])/norm2(a[[q]])
        Y[, q] <- pm(A[[q]], a[[q]], na.rm = na.rm)
      }

    # check for convergence of the SGCCA algorithm
    ifelse((mode(scheme) != "function"),
           {g <- function(x) switch(scheme,
                                    horst = x,
                                    factorial = x**2,
                                    centroid = abs(x))
            crit[iter] <- sum(C*g(cov2(Y, bias = bias)))
            },
           crit[iter] <- sum(C*scheme(cov2(Y, bias = bias)))
    )

    # Print out intermediate fit

    if (verbose & (iter %% 1)==0)
     cat(" Iter: ", formatC(iter,width=3, format="d"),
         " Fit: ", formatC(crit[iter], digits=8, width=10, format="f"),
         " Dif: ", formatC(crit[iter]-crit_old, digits=8, width=10, format="f"),
         "\n")

    stopping_criteria = c(drop(crossprod(Reduce("c", mapply("-", a, a_old))))
                          , abs(crit[iter]-crit_old))

    if ( any(stopping_criteria < tol) | (iter > 1000))
      break

    crit_old = crit[iter]
    a_old <- a
    iter <- iter + 1
  }


  if (iter > 1000) warning("The SGCCA algorithm did not converge after 1000
                           iterations.")
  if(iter<1000 & verbose) cat("The SGCCA algorithm converged to a stationary
                              point after", iter-1, "iterations \n")
  if (verbose) plot(crit, xlab = "iteration", ylab = "criteria")

  for (q in 1:J) if(sum(a[[q]]!=0) <= 1)
      {
        if(!quiet)
        {
            warning(sprintf("Deflation failed because only one variable was
                            selected for block #",q))
        }
     }

  AVE_inner  <- sum(C*cor(Y)^2/2)/(sum(C)/2) # AVE inner model

  result <- list(Y = Y, a = a, crit = crit[which(crit != 0)],
                 AVE_inner = AVE_inner)
  return(result)
}

