#' SGCCA extends RGCCA to address the issue of variable selection. Specifically, 
#' RGCCA is combined with an L1-penalty that gives rise to Sparse GCCA (SGCCA) 
#' which is implemented in the function sgcca().
#' Given \eqn{J} matrices \eqn{X_1, X_2, ..., X_J}, that represent 
#' \eqn{J} sets of variables observed on the same set of \eqn{n} individuals. 
#' The matrices \eqn{X_1, X_2, ..., X_J} must have the same number of rows, but 
#' may (and usually will) have different numbers of columns. Blocks are not 
#' necessarily fully connected within the SGCCA framework. Hence the use of 
#' SGCCA requires the construction (user specified) of a design matrix (\eqn{C}) 
#' that characterizes the connections between blocks. Elements of the symmetric 
#' design matrix \eqn{C = (c_{jk})} are equal to 1 if block \eqn{j} and block 
#' \eqn{k} are connected, and 0 otherwise. The SGCCA algorithm is very similar 
#' to the RGCCA algorithm and keeps the same monotone convergence properties 
#' (i.e. the bounded criteria to be maximized increases at each step of the 
#' iterative procedure and hits at convergence a stationary point).
#' Moreover, using a deflation strategy, sgcca() enables the computation of 
#' several SGCCA block components (specified by ncomp) for each block. Block 
#' components for each block are guaranteed to be orthogonal when using this 
#' deflation strategy. The so-called symmetric deflation is considered in this 
#' implementation, i.e. each block is deflated with respect to its own 
#' component. Moreover, we stress that the numbers of components per block 
#' could differ from one block to another. 
#' @inheritParams select_analysis
#' @inheritParams rgccaNa
#' @inheritParams rgccad
#' @param sparsity Either a \eqn{1*J} vector or a \eqn{max(ncomp) * J} matrix 
#' encoding the L1 constraints applied to the outer weight vectors. The amount 
#' of sparsity varies between \eqn{1/sqrt(p_j)} and 1 (larger values of sparsity 
#' correspond to less penalization). If sparsity is a vector, L1-penalties are 
#' the same for all the weights corresponding to the same block but different 
#' components: 
#' \deqn{for all h, |a_{j,h}|_{L_1} \le c_1[j] \sqrt{p_j},}
#' with \eqn{p_j} the number of variables of \eqn{X_j}.
#' If sparsity is a matrix, each row \eqn{h} defines the constraints applied to 
#' the weights corresponding to components \eqn{h}:
#' \deqn{for all h, |a_{j,h}|_{L_1} \le c_1[h,j] \sqrt{p_j}.} It can be 
#' estimated by using \link{rgcca_permutation}.
#' @return \item{Y}{A list of \eqn{J} elements. Each element of \eqn{Y} is a 
#' matrix that contains the analysis components for the corresponding block.}
#' @return \item{a}{A list of \eqn{J} elements. Each element of \eqn{a} is a 
#' matrix that contains the outer weight vectors for each block.}
#' @return \item{astar}{A list of \eqn{J} elements. Each element of astar is a 
#' matrix defined as Y[[j]][, h] = A[[j]]\%*\%astar[[j]][, h]}
#' @return \item{crit}{A vector of integer that contains for each component the 
#' values of the analysis criteria across iterations.}
#' @return \item{AVE}{A list of numerical values giving the indicators of model 
#' quality based on the Average Variance Explained (AVE): AVE(for each block), 
#' AVE(outer model), AVE(inner model).}
#' @references Tenenhaus, A., Philippe, C., Guillemot, V., Le Cao, K. A., 
#' Grill, J., and Frouin, V. , "Variable selection for generalized canonical 
#' correlation analysis.," Biostatistics, vol. 15, no. 3, pp. 569-583, 2014. 
#' @title Variable Selection For Generalized Canonical Correlation Analysis 
#' (SGCCA)
#' @examples
#'
#' #############
#' # Example 1 #
#' #############
#' \dontrun{
#' # Download the dataset's package at http://biodev.cea.fr/sgcca/.
#' # --> gliomaData_0.4.tar.gz
#' 
#' require(gliomaData)
#' data(ge_cgh_locIGR)
#' 
#' A <- ge_cgh_locIGR$multiblocks
#' Loc <- factor(ge_cgh_locIGR$y)
#' levels(Loc) <- colnames(ge_cgh_locIGR$multiblocks$y)
#' C <-  matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3)
#' tau = c(1, 1, 0)
#' 
#' # rgcca algorithm using the dual formulation for X1 and X2 
#' # and the dual formulation for X3
#' A[[3]] = A[[3]][, -3]
#' result.rgcca = rgccad(A, C, tau, ncomp = c(2, 2, 1), 
#' scheme = "factorial", verbose = TRUE)
#' # sgcca algorithm
#' result.sgcca = sgcca(A, C, sparsity = c(.071,.2, 1), ncomp = c(2, 2, 1),
#'                      scheme = "centroid", verbose = TRUE)
#' 
#' ############################
#' # plot(y1, y2) for (RGCCA) #
#' ############################
#' layout(t(1:2))
#' plot(result.rgcca$Y[[1]][, 1], result.rgcca$Y[[2]][, 1], 
#'      col = "white", xlab = "Y1 (GE)", ylab = "Y2 (CGH)", 
#'      main = "Factorial plan of RGCCA")
#' text(result.rgcca$Y[[1]][, 1], result.rgcca$Y[[2]][, 1], 
#'      Loc, col = as.numeric(Loc), cex = .6)
#' plot(result.rgcca$Y[[1]][, 1], result.rgcca$Y[[1]][, 2], 
#'      col = "white", xlab = "Y1 (GE)", 
#'      ylab = "Y2 (GE)", main = "Factorial plan of RGCCA")
#' text(result.rgcca$Y[[1]][, 1], result.rgcca$Y[[1]][, 2], 
#'      Loc, col = as.numeric(Loc), cex = .6)
#' 
#' ############################
#' # plot(y1, y2) for (SGCCA) #
#' ############################
#' layout(t(1:2))
#' plot(result.sgcca$Y[[1]][, 1], result.sgcca$Y[[2]][, 1], 
#'      col = "white", xlab = "Y1 (GE)", 
#'      ylab = "Y2 (CGH)", main = "Factorial plan of SGCCA")
#' text(result.sgcca$Y[[1]][, 1], result.sgcca$Y[[2]][, 1], 
#'      Loc, col = as.numeric(Loc), cex = .6)
#' 
#' plot(result.sgcca$Y[[1]][, 1], result.sgcca$Y[[1]][, 2], 
#'      col = "white", xlab = "Y1 (GE)", 
#'      ylab = "Y2 (GE)", main = "Factorial plan of SGCCA")
#' text(result.sgcca$Y[[1]][, 1], result.sgcca$Y[[1]][, 2], 
#'      Loc, col = as.numeric(Loc), cex = .6)
#' 
#' # sgcca algorithm with multiple components and different 
#' L1 penalties for each components 
#' # (-> sparsity is a matrix)
#' init = "random"
#' result.sgcca = sgcca(A, C, 
#'                      sparsity = matrix(c(.071,.2, 1, 0.06, 0.15, 1), 
#'                      nrow = 2, byrow = TRUE),
#'                      ncomp = c(2, 2, 1), scheme = "factorial", 
#'                      scale = TRUE, bias = TRUE, 
#'                      init = init, verbose = TRUE)
#' # number of non zero elements per dimension
#' apply(result.sgcca$a[[1]], 2, function(x) sum(x!=0)) 
#'      #(-> 145 non zero elements for a11 and 107 non zero elements for a12)
#' apply(result.sgcca$a[[2]], 2, function(x) sum(x!=0)) 
#'      #(-> 85 non zero elements for a21 and 52 non zero elements for a22)
#' init = "svd"
#' result.sgcca = sgcca(A, C, sparsity = matrix(c(.071,.2, 1, 0.06, 0.15, 1), 
#'                      nrow = 2, byrow = TRUE),
#'                      ncomp = c(2, 2, 1), scheme = "factorial", scale = TRUE, 
#'                      bias = TRUE, 
#'                      init = init, verbose = TRUE)
#' }
#'@export sgcca

sgcca <- function (A, C = 1-diag(length(A)), sparsity = rep(1, length(A)), 
                   ncomp = rep(1, length(A)), scheme = "centroid", scale = TRUE, 
                   init = "svd", bias = TRUE, tol = .Machine$double.eps, 
                   verbose = FALSE, scale_block = TRUE, prescaling = FALSE,
                   quiet = FALSE){
  
  ndefl <- ncomp-1
  N <- max(ndefl)
  J <- length(A)
  pjs <- sapply(A,NCOL)
  nb_ind <- NROW(A[[1]])
  AVE_X = list()
  AVE_outer <- rep(NA,max(ncomp))
  
  if ( any(ncomp < 1) ) stop_rgcca("One must compute at least one component per 
                                   block!")
  if (any(ncomp-pjs > 0)) stop_rgcca("For each block, choose a number of 
                                     components smaller than the number of 
                                     variables!")
  
  if (is.vector(sparsity)){
    if (any(sparsity < 1/sqrt(pjs) | sparsity > 1 ))
      stop_rgcca("L1 constraints (sparsity) must vary between 1/sqrt(p_j) and 1.")
  }
  
  if (is.matrix(sparsity)){
    if (any(apply(sparsity, 1, function(x) any(x < 1/sqrt(pjs)))))
      stop_rgcca("L1 constraints (sparsity) must vary between 1/sqrt(p_j) 
                 and 1.")
  }
  
###################################################
 
  if (mode(scheme) != "function") {
    if ((scheme != "horst" ) & (scheme != "factorial") & (scheme != "centroid")) {
      stop_rgcca("Choose one of the three following schemes: horst, centroid, 
                 factorial or design the g function")
    }
    if (verbose) cat("Computation of the SGCCA block components based on the", 
                     scheme, "scheme \n")
  }
  if (mode(scheme) == "function" & verbose) {
    cat("Computation of the SGCCA block components based on the g scheme \n")
  }
  
  
    #-------------------------------------------------------
    if(!prescaling)
        A=scaling(A, scale = scale, bias = bias, scale_block = scale_block)

    ####################################
    # sgcca with 1 component per block #
    ####################################
          
    # ndefl number of deflation per block
    ndefl <- ncomp-1
    N <- max(ndefl)
    J <- length(A)
    pjs <- sapply(A,NCOL)
    nb_ind <- NROW(A[[1]])
    AVE_X = list()
    AVE_outer <- rep(NA,max(ncomp))
    if (N == 0) {
        result <- sgccak(A, C, sparsity, scheme, init = init, bias = bias, 
                         tol = tol, verbose = verbose, quiet = quiet)
        # No deflation (No residual matrices generated).
        Y <- NULL
        for (b in 1:J) Y[[b]] <- result$Y[,b, drop = FALSE]
        #Average Variance Explained (AVE) per block
        for (j in 1:J) AVE_X[[j]] =  mean(cor(A[[j]], Y[[j]],
                                              use="pairwise.complete.obs")^2, 
                                          na.rm=TRUE)
        
        #AVE outer 
        AVE_outer <- sum(pjs * unlist(AVE_X))/sum(pjs)
        
        AVE <- list(AVE_X = AVE_X, 
                    AVE_outer = AVE_outer,
                    AVE_inner = result$AVE_inner)
        
        a <- lapply(result$a, cbind)
        
        for (b in 1:J) {
          rownames(a[[b]]) = colnames(A[[b]])  
          rownames(Y[[b]]) = rownames(A[[b]])
          colnames(Y[[b]]) = "comp1"
        }
        
        out <- list(Y=Y, a=a, astar=a, 
                    C=C, scheme=scheme, sparsity=sparsity, ncomp=ncomp,
                    crit = result$crit[length(result$crit)],
                    AVE = AVE, blocks = A)
        class(out) <- "sgcca"
        return(out)
    }
    
    ##################
    # Initialization #
    ##################
    
    Y <- NULL
    R <- A
    P <- a <- astar <- NULL
    crit <- list()
    AVE_inner <- rep(NA,max(ncomp))
    
    for (b in 1:J) P[[b]] <- a[[b]] <- astar[[b]] <- matrix(NA,pjs[[b]],N+1)
    for (b in 1:J) Y[[b]] <- matrix(NA,nb_ind,N+1)
    
    ##############################################
    #               If any ncomp > 1             #
    #      Determination of SGCCA components     #
    ##############################################
    

    for (n in 1:N) {
      if (verbose) cat(paste0("Computation of the SGCCA block components #", n, 
                              " is under progress... \n"))
      if(is.vector(sparsity)){
        sgcca.result <- sgccak(R, C, sparsity = sparsity , scheme=scheme, 
                               init = init, bias = bias, tol = tol, 
                               verbose = verbose, quiet = quiet)
      } else{
        sgcca.result <- sgccak(R, C, sparsity = sparsity[n, ] , scheme = scheme, 
                               init = init, bias = bias, tol = tol,  
                               verbose = verbose, quiet = quiet)
      }
      AVE_inner[n] <- sgcca.result$AVE_inner
      crit[[n]] <- sgcca.result$crit
      
      
      for (b in 1:J) Y[[b]][,n] <- sgcca.result$Y[ ,b]
      for (q in which(n <ndefl)) if(sum(sgcca.result$a[[q]]!=0) <= 1)
     { 
        if(!quiet)
        {
            warning(sprintf("Deflation failed because only one variable was 
                            selected for block ",q,"! \n"))
            
        }
     }
	    defla.result <- defl.select(sgcca.result$Y, R, ndefl, n, nbloc = J)
      R <- defla.result$resdefl
      for (b in 1:J) {
        P[[b]][,n] <- defla.result$pdefl[[b]]
        a[[b]][,n] <- sgcca.result$a[[b]]
      }
      
      if (n==1) {
        for (b in 1:J) astar[[b]][,n] <- sgcca.result$a[[b]]
      } 
      else {
        for (b in 1:J) astar[[b]][,n] <- sgcca.result$a[[b]] -  
            astar[[b]][,(1:n-1),drop=F] %*% 
            drop( t(a[[b]][,n]) %*% P[[b]][,1:(n-1),drop=F] ) 
      }     
    }
    if (verbose) cat(paste0("Computation of the SGCCA block components #", N+1, " 
                            is under progress...\n"))
    if(is.vector(sparsity)) {
      sgcca.result <- sgccak(R, C, sparsity = sparsity, scheme=scheme, 
                             init = init, bias = bias, tol = tol, 
                             verbose = verbose, quiet = quiet)
    } else{
      sgcca.result <- sgccak(R, C, sparsity = sparsity[N+1, ], scheme=scheme, 
                             init = init, bias = bias, tol = tol, 
                             verbose = verbose, quiet = quiet)
    }
    AVE_inner[max(ncomp)] <- sgcca.result$AVE_inner
    
    crit[[N+1]] <- sgcca.result$crit
    for (b in 1:J) {
      Y[[b]][,N+1]     <- sgcca.result$Y[, b]
      a[[b]][,N+1]     <- sgcca.result$a[[b]]
      astar[[b]][,N+1] <- sgcca.result$a[[b]] -  astar[[b]][,(1:N),drop=F] %*% 
                              drop(t(a[[b]][, N+1]) %*% P[[b]][,1:N,drop=F]) 
      rownames(a[[b]]) = rownames(astar[[b]]) = colnames(A[[b]])  
      rownames(Y[[b]]) = rownames(A[[b]])
      colnames(Y[[b]]) = paste0("comp", 1:max(ncomp))
    }
    
    #Average Variance Explained (AVE) per block
    for (j in 1:J) AVE_X[[j]] =  apply(cor(A[[j]], Y[[j]],
                                           use="pairwise.complete.obs")^2, 2, 
                                       function(x) {return(mean(x,is.na=TRUE))})

    #AVE outer 
    outer = matrix(unlist(AVE_X), nrow = max(ncomp))
    for (j in 1:max(ncomp)) 
      AVE_outer[j] <- sum(pjs * outer[j, ], na.rm=T)/sum(pjs)
    
    Y = shave.matlist(Y, ncomp)
    AVE_X = shave.veclist(AVE_X, ncomp)
    
    AVE <- list(AVE_X = AVE_X, 
                AVE_outer = AVE_outer,
                AVE_inner = AVE_inner)
    
  out <- list(Y = shave.matlist(Y, ncomp),
                a = shave.matlist(a, ncomp), 
                astar = shave.matlist(astar, ncomp),
                crit = crit,
                AVE = AVE, blocks=A
                )

    class(out) <- "sgcca"
    return(out)
    
}
