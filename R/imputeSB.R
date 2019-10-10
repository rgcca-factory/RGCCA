#'Impute with superblock method
#'
#'This method is used for the implementation of EM algorithm for missing data
#'
#' @param A A list of J blocks
#' @param tau A vector of tau values with the same length as A
#' @param tol The stopping value for convergence.
#' @param graph if graph = TRUE, 
#' @param ncomp vector containing the number of components per block in RGCCA
#' @param naxis number of component to select in the superblock for the estimation of missing data
#' @param  scale  If scale = TRUE, each block is standardized to zero means and unit variances (default: TRUE).
#' @param sameBlockWeight A logical value indicating if the different blocks should have the same weight in the analysis (default, sameBlockWeight=TRUE)
#' @param bias A logical value indicating if variance should be biased or not
#' @return \item{A} A list of blocks imputed
#' @return \item{crit} Convergence criterion : abs(1-obj_k/obj_{k-1})
#' @return \item{obj} Vector containing the mean square error between the predict values and the original non missing values at each iteration
#' @title imputeSB: impute with superblock method
#' @examples 
#'  data();...

# TODO: tau did not have a default value
# TODO: scheme in par
imputeSB <- function(
    A,
    tau = rep(1, length(A)),
    ni = 50,
    tol = 1e-08,
    graph = FALSE,
    ncomp = NULL,
    naxis = 1,
    scale = TRUE,
    sameBlockWeight = TRUE,
    bias = TRUE) {
    
    # listWithoutNA
    nvar <- sapply(A, NCOL)
    indNA <- which(is.na(do.call(cbind, A)), arr.ind = TRUE)

    # to get the superblock X1NA, with missing values imputed by the colmeans
    Alist <- lapply(A, function(X) {
        if (is.list(X))
            X <- as.matrix(X)

        if (is.matrix(X)) {
            m <- apply(X, 2, function(x) mean(x, na.rm = TRUE))
            indNA <- which(is.na(X), arr.ind = TRUE)
            X[indNA] <- m[indNA[, 2]]
            return(X)
        }
        if (is.vector(X)) {
            m <- mean(X, na.rm = TRUE)
            X[is.na(X)] <- m
            return(X)
        }
    })

    X1NA <- Reduce(cbind, Alist)
    # center and normalized ???
    D <- matrix(1, dim(X1NA)[1], dim(X1NA)[2])
    X2NA <- scale2(X1NA, scale = scale, bias = TRUE)

    if (sameBlockWeight) {
        group <- unlist(lapply(A, "NCOL"))
        debutBlock <- c(1, 1 + cumsum(group)[1:(length(group) - 1)])
        finBlock <- cumsum(group)
        for (u in 1:length(finBlock)) {
            var_group <- sum(apply(X2NA[, debutBlock[u]:finBlock[u]], 2, "cov2"))
            D[, debutBlock[u]:finBlock[u]] <- 1/sqrt(var_group)
        }
    }
    # initialization
    i <- 1
    diff <- objective <- old <- criterion <- list()
    criterionFinal <- objectiveFinal <- list()
    continue <- TRUE
    
    J <- length(A)
    C2 <- matrix(1, J + 1, J + 1)
    C2[1:J, 1:J] <- 0
    C2[J + 1, J + 1] <- 0
    tau2 <- c(tau, 0)  
    # mode A for the blocks and mode B for superblock (design for CPCA2 )

    if (is.null(ncomp))
        ncomp2 <- c(rep(1, J), naxis)
    else
        ncomp2 <- c(ncomp, naxis)

    critRGCCA <- c()
    old[[1]] <- -Inf
    
# Reconstruction ! ! !
    Xhat = (y%*%t(a))*stdev*(1/D) + moy
    
    # 
    X1NA[indNA] = Xhat[indNA]
    Alist=list()
    for(j in 1:length(nvar))
    {
      if(j==1){sel=1:nvar[1]}else{debut=sum(nvar[1:(j-1)])+1;fin=debut+(nvar[j]-1);sel=debut:fin}
      Alist[[j]]=as.matrix(X1NA[,sel])
      colnames( Alist[[j]])=colnames(A[[j]])
    }
    
    if (graph) {
        vec <- rep(NA, i - 1)
        for (k in 1:(i - 1)) 
            vec[k] <- criterion[[k]]

        # x11()
        png("objective.png")
        plot(vec[-1], xlab = "iteration", ylab = "objective", pch = 16)
        dev.off()
        
        png("criterion.png")
        plot(critRGCCA, xlab = "iteration", ylab = "criterion", pch = 16)
        dev.off()
    }
    if (i > ni)
    {
      continue <- FALSE
      warning(paste("The RGCCA imputation algorithm did not converge after ", ni, " iterations"))
    }
    i<- i + 1
  }	
  
   if(graph)
   {   vec=rep(NA,i-1)
   for(k in 1:(i-1)){vec[k]=criterion[[k]]}
  # x11()
   png("objective.png")
   plot(vec[-1], xlab = "iteration", ylab = "objective",pch=16)
   dev.off()
  
   png("criterion.png")
   plot(critRGCCA, xlab = "iteration", ylab = "criterion",pch=16)
   dev.off()
  }

  return(list(A=Alist,crit=unlist(critRGCCA),stab=unlist(criterion)))
}
