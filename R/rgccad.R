#' Regularized Generalized Canonical Correlation Analysis (RGCCA) is a generalization
#' of regularized canonical correlation analysis to three or more sets of variables. 
#' Given \eqn{J} matrices \eqn{\mathbf{X_1}, \mathbf{X_2}, ..., \mathbf{X_J}} that represent 
#' \eqn{J} sets of variables observed on the same set of \eqn{n} individuals. The matrices 
#' \eqn{\mathbf{X_1}, \mathbf{X_2}, ..., \mathbf{X_J}} must have the same number of rows, 
#' but may (and usually will) have different numbers of columns. The aim of RGCCA is to study 
#' the relationships between these \eqn{J} blocks of variables. It constitutes a general 
#' framework for many multi-block data analysis methods. It combines the power of 
#' multi-block data analysis methods (maximization of well identified criteria) 
#' and the flexibility of PLS path modeling (the researcher decides which blocks 
#' are connected and which are not). Hence, the use of RGCCA requires the construction 
#' (user specified) of a design matrix, (\eqn{\mathbf{C}}), that characterize 
#' the connections between blocks. Elements of the (symmetric) design matrix \eqn{\mathbf{C} = (c_{jk})} 
#' is equal to 1 if block \eqn{j} and block \eqn{k} are connected, and 0 otherwise.
#' The objective is to find a fixed point of the stationary equations related to the RGCCA optimization 
#' problem. The function rgccad() implements a monotonically convergent algorithm (i.e. the bounded
#' criteria to be maximized increases at each step of the iterative procedure) that is very 
#' similar to the PLS algorithm proposed by Herman Wold. Moreover, depending on the 
#' dimensionality of each block \eqn{\mathbf{X}_j}, \eqn{j = 1, \ldots, J}, the primal (when \eqn{n > p_j}) algorithm or 
#' the dual (when \eqn{n < p_j}) algorithm is used (see Tenenhaus et al. 2013). 
#' Moreover, by deflation strategy, rgccad() allow to compute several RGCCA block
#' components (specified by ncomp) for each block. Block components of each block are guaranteed to 
#' be orthogonal with the use of the deflation. The so-called symmetric deflation is considered in
#' this implementation, i.e. each block is deflated with respect to its own component.
#' It should be noted that the numbers of components per block can differ from one block to another. 
#' @param A  A list that contains the \eqn{J} blocks of variables \eqn{\mathbf{X_1}, \mathbf{X_2}, ..., \mathbf{X_J}}.
#' @param C  A design matrix that describes the relationships between blocks (default: complete design).
#' @param tau tau is either a \eqn{1 \times J} vector or a \eqn{\mathrm{max}(ncomp) \times J} matrix, and contains the values 
#' of the shrinkage parameters (default: tau = 1, for each block and each dimension).
#' If tau = "optimal" the shrinkage paramaters are estimated for each block and each dimension using the Schafer and Strimmer (2005)
#' analytical formula . If tau is a \eqn{1\times J} numeric vector, tau[j] is identical across the dimensions of block \eqn{\mathbf{X}_j}. 
#' If tau is a matrix, tau[k, j] is associated with \eqn{\mathbf{X}_{jk}} (\eqn{k}th residual matrix for block \eqn{j})
#' @param ncomp  A \eqn{1 \times J} vector that contains the numbers of components for each block (default: rep(1, length(A)), which gives one component per block.)
#' @param scheme The value is "horst", "factorial", "centroid" or the g function (default: "centroid").
#' @param scale  If scale = TRUE, each block is standardized to zero means and unit variances (default: TRUE).
#' @param verbose  If verbose = TRUE, the progress will be report while computing (default: TRUE).
#' @param init The mode of initialization to use in RGCCA algorithm. The alternatives are either by Singular Value Decompostion ("svd") or random ("random") (Default: "svd").
#' @param bias A logical value for biaised or unbiaised estimator of the var/cov (default: bias = TRUE).
#' @param tol The stopping value for convergence.
#' @param sameBlockWeight A logical value indicating if the different blocks should have the same weight in the analysis (default, sameBlockWeight=TRUE)
#' @param na.rm If TRUE, runs rgcca only on available data.
#' @param estimateNA If TRUE, missing values are estimated during the RGCCA calculation
#' @return \item{Y}{A list of \eqn{J} elements. Each element of \eqn{Y} is a matrix that contains the RGCCA components for the corresponding block.}
#' @return \item{a}{A list of \eqn{J} elements. Each element of \eqn{a} is a matrix that contains the outer weight vectors for each block.}
#' @return \item{astar}{A list of \eqn{J} elements. Each element of astar is a matrix defined as Y[[j]][, h] = A[[j]]\%*\%astar[[j]][, h].}
#' @return \item{tau}{A vector or matrix that contains the values of the shrinkage parameters applied to each block and each dimension (user specified).}
#' @return \item{call}{Call of the function}#' @return \item{crit}{A vector that contains the values of the objective function at each iterations.}
#' @return \item{crit}{A vector that contains the values of the criteria across iterations.}
#' @return \item{mode}{A \eqn{1 \times J} vector that contains the formulation ("primal" or "dual") applied to each of the \eqn{J} blocks within the RGCCA alogrithm} 
#' @return \item{AVE}{indicators of model quality based on the Average Variance Explained (AVE): AVE(for one block), AVE(outer model), AVE(inner model).}
#' @references Tenenhaus A. and Tenenhaus M., (2011), Regularized Generalized Canonical Correlation Analysis, Psychometrika, Vol. 76, Nr 2, pp 257-284.
#' @references Tenenhaus A. et al., (2013), Kernel Generalized Canonical Correlation Analysis, submitted.
#' @references Schafer J. and Strimmer K., (2005), A shrinkage approach to large-scale covariance matrix estimation and implications for functional genomics. Statist. Appl. Genet. Mol. Biol. 4:32.
#' @title Regularized Generalized Canonical Correlation Analysis (RGCCA) 
#' @examples
#' #############
#' # Example 1 #
#' #############
#' data(Russett)
#' X_agric =as.matrix(Russett[,c("gini","farm","rent")])
#' X_ind = as.matrix(Russett[,c("gnpr","labo")])
#' X_polit = as.matrix(Russett[ , c("demostab", "dictator")])
#' A = list(X_agric, X_ind, X_polit)
#' #Define the design matrix (output = C) 
#' C = matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3)
#' result.rgcca = rgccad(A, C, tau = c(1, 1, 1), scheme = "factorial", scale = TRUE)
#' lab = as.vector(apply(Russett[, 9:11], 1, which.max))
#' plot(result.rgcca$Y[[1]], result.rgcca$Y[[2]], col = "white", 
#'      xlab = "Y1 (Agric. inequality)", ylab = "Y2 (Industrial Development)")
#' text(result.rgcca$Y[[1]], result.rgcca$Y[[2]], Russett[, 1], col = lab, cex = .7)
#'
#' ############################################
#' # Example 2: RGCCA and mutliple components #
#' ############################################
#' ############################
#' # plot(y1, y2) for (RGCCA) #
#' ############################
#' result.rgcca = rgccad(A, C, tau = rep(1, 3), ncomp = c(2, 2, 1),
#'                      scheme = "factorial", verbose = TRUE)
# 'layout(t(1:2))
#' plot(result.rgcca$Y[[1]][, 1], result.rgcca$Y[[2]][, 1], col = "white", xlab = "Y1 (GE)", 
#'     ylab = "Y2 (CGH)", main = "Factorial plan of RGCCA")
#' text(result.rgcca$Y[[1]][, 1], result.rgcca$Y[[2]][, 1], Russett[, 1], col = lab, cex = .6)
#' plot(result.rgcca$Y[[1]][, 1], result.rgcca$Y[[1]][, 2], col = "white", xlab = "Y1 (GE)", 
#'     ylab = "Y2 (GE)", main = "Factorial plan of RGCCA")
#' text(result.rgcca$Y[[1]][, 1], result.rgcca$Y[[1]][, 2], Russett[, 1], col = lab, cex = .6)
#'
#' ######################################
#' # example 3: RGCCA and leave one out #
#' ######################################
#' Ytest = matrix(0, 47, 3)
#' result.rgcca = rgccad(A, C, tau = rep(1, 3), ncomp = rep(1, 3),
#'                      scheme = "factorial", verbose = TRUE)
#'                      
#' for (i in 1:nrow(Russett)){
#'  B = lapply(A, function(x) x[-i, ])
#'  B = lapply(B, scale2)
#'  resB = rgccad(B, C, tau = rep(1, 3), scheme = "factorial", scale = FALSE, verbose = FALSE)
#'  #  look for potential conflicting sign among components within the loo loop.
#'  for (k in 1:length(B)){
#'    if (cor(result.rgcca$a[[k]], resB$a[[k]]) >= 0) 
#'      resB$a[[k]] = resB$a[[k]] else resB$a[[k]] = -resB$a[[k]]
#'  }
#'  Btest =lapply(A, function(x) x[i, ])
#'  Btest[[1]]=(Btest[[1]]-attr(B[[1]],"scaled:center"))/(attr(B[[1]],"scaled:scale"))
#'  Btest[[2]]=(Btest[[2]]-attr(B[[2]],"scaled:center"))/(attr(B[[2]],"scaled:scale"))
#'  Btest[[3]]=(Btest[[3]]-attr(B[[3]],"scaled:center"))/(attr(B[[3]],"scaled:scale"))
#'  Ytest[i, 1] = Btest[[1]]%*%resB$a[[1]]
#'  Ytest[i, 2] = Btest[[2]]%*%resB$a[[2]]
#'  Ytest[i, 3] = Btest[[3]]%*%resB$a[[3]]
#' }
#' lab = apply(Russett[, 9:11], 1, which.max)
#' plot(result.rgcca$Y[[1]], result.rgcca$Y[[2]], col = "white", 
#'      xlab = "Y1 (Agric. inequality)", ylab = "Y2 (Ind. Development)")
#' text(result.rgcca$Y[[1]], result.rgcca$Y[[2]], Russett[, 1], col = lab)
#' text(Ytest[, 1], Ytest[, 2], substr(Russett[, 1], 1, 1), col = lab)
#' @export rgccad
#' @importFrom grDevices dev.off png rainbow
#' @importFrom graphics abline axis close.screen grid legend lines par points rect screen segments split.screen text
#' @importFrom stats binomial glm lm predict sd var weighted.mean
#' @importFrom utils read.table write.table

#' @importFrom grDevices graphics.off

rgccad=function (A, C = 1 - diag(length(A)), tau = rep(1, length(A)),  ncomp = rep(1, length(A)), scheme = "centroid", scale = TRUE,   init = "svd", bias = TRUE, tol = 1e-08, verbose = TRUE,sameBlockWeight=TRUE,na.rm=TRUE,estimateNA="no")
{
  shave.matlist <- function(mat_list, nb_cols) mapply(function(m,nbcomp) m[, 1:nbcomp, drop = FALSE], mat_list, nb_cols, SIMPLIFY = FALSE)
  shave.veclist <- function(vec_list, nb_elts) mapply(function(m, nbcomp) m[1:nbcomp], vec_list, nb_elts, SIMPLIFY = FALSE)
  A0=A
  #  call = match.call()
  call=list(A=A, C = C, tau = tau,  ncomp = ncomp, scheme = scheme, scale = scale,   init = init, bias = bias, tol =tol, verbose = verbose,sameBlockWeight=sameBlockWeight,na.rm=na.rm,estimateNA=estimateNA)
  if (any(ncomp < 1)) {stop("Compute at least one component per block!")}	
  pjs <- sapply(A, NCOL) #nombre de variables par bloc
  varij <- sapply(A,function(x){covarMat=cov2(x,bias=bias);varianceBloc=sum(diag(covarMat));return(varianceBloc)})
 # print("varij")
 # print(varij)
  nb_row <- NROW(A[[1]]) #nombre de lignes
  
  # Verifications des commandes entrees par l'utilisateur
  if (any(ncomp - pjs > 0))  # le nombre de composantes doit etre inferieur au nombre de variables
    stop("For each block, choose a number of components smaller than the number of variables!")
  if (mode(scheme) != "function") {
    if ((scheme != "horst") & (scheme != "factorial") & (scheme != 
                                                         "centroid")) {
      stop("Choose one of the three following schemes: horst, centroid, factorial or design the g function")
    }
    if (verbose) 
      cat("Computation of the RGCCA block components based on the", 
          scheme, "scheme \n")
  }
  if (mode(scheme) == "function" & verbose) {
    cat("Computation of the RGCCA block components based on the g scheme \n")
  }
  if (!is.numeric(tau) & verbose) { #tau peut etre estime
    cat("Optimal Shrinkage intensity paramaters are estimated \n")
  }
  else {
    if (is.numeric(tau) & verbose) {
      cat("Shrinkage intensity paramaters are chosen manually \n")
    }
  }

      if (scale == TRUE) 
      {
          
          A = lapply(A, function(x) scale2(x,scale=TRUE, bias = bias)) # le biais indique si on recherche la variance biaisee ou non
          if(sameBlockWeight)
          {
              A = lapply(A, function(x) 
                  {
                    y=x/sqrt(NCOL(x));
                    return(y)
                  } 
            )
          }
          # on divise chaque bloc par la racine du nombre de variables pour avoir chaque poids pour le meme bloc
      }
      if (scale == FALSE)
      { 
          
          A = lapply(A, function(x) scale2(x, scale=FALSE, bias = bias)) 
          if(sameBlockWeight)
          {
              A = lapply(A, function(x) 
                    {
                        covarMat=cov2(x,bias=bias);
                        varianceBloc=sum(diag(covarMat)); 
                        return(x/sqrt(varianceBloc))
                    })
          }
          
      }
      
#  }

  # Superblock option
  if(!is.matrix(C)&& C=="superblock")
  {
    A=c(A,list(do.call(cbind,A)))
    C=matrix(0,length(A),length(A))
    C[length(A),1:(length(A)-1)]=1
    C[1:(length(A)-1),length(A)]=1
    tau=c(tau,1)
    pjs=c(pjs,sum(pjs))
    ncomp=c(ncomp,ncomp[1])
    #tau=c(tau,0)
  }
  AVE_X = list() 
  AVE_outer <- vector()
  ndefl <- ncomp - 1
  N <- max(ndefl)
  nb_ind <- NROW(A[[1]])
  J <- length(A)
  # si le nombre d'individu est inferieur au nombre de variables: primal, sinon dual
  primal_dual = rep("primal", J)
  primal_dual[which(nb_row < pjs)] = "dual"
  # cas ou le nombre de composantes
  if (N == 0) 
  { # cas ou on n'a qu'un axe a calculer par bloc
    
    result <- rgccak(A, C, tau = tau, scheme = scheme, init = init, bias = bias, tol = tol, verbose = verbose,na.rm=na.rm,estimateNA=estimateNA,sameBlockWeight=sameBlockWeight,scale=scale)
    if(estimateNA%in%c("iterative","first","lebrusquet","superblock"))
    {
      A<-result$call$A
    }
    # result contient le resultat de rgcca
    Y <- NULL 
    for (b in 1:J) Y[[b]] <- result$Y[, b, drop = FALSE]
  
        for (j in 1:J)
        {
           # AVE_X[[j]] = mean(cor(A[[j]], Y[[j]],use="pairwise.complete.obs")^2,na.rm=TRUE)
            AVE_X[[j]]=diag(cov(Y[[j]]))/sum(diag(cov(A[[j]] )))
        } 
        AVE_outer <- sum(pjs * unlist(AVE_X))/sum(pjs) 
   
    AVE <- list(AVE_X = AVE_X, AVE_outer = AVE_outer, AVE_inner = result$AVE_inner)
    a <- lapply(result$a, cbind)
    for (b in 1:J) {
      rownames(a[[b]]) = colnames(A[[b]])
      rownames(Y[[b]]) = rownames(A[[b]])
      colnames(Y[[b]]) = "comp1"
    }
    out <- list(Y = Y, a = a, astar = a, C = C, tau = result$call$tau,  scheme = scheme, ncomp = ncomp, crit = result$crit, primal_dual = primal_dual, AVE = AVE,A=A0,call=call)
    if(estimateNA %in% c("iterative","first","superblock","lebrusquet")){out[["imputedA"]]=A}
    
    class(out) <- "rgccad"
    
    return(out)
  }
  Y <- NULL
  crit = list()
  AVE_inner <- rep(NA, max(ncomp))
  R <- A
  P <- a <- astar <- NULL
  if (is.numeric(tau)) 
    tau_mat = tau
  else tau_mat = matrix(NA, max(ncomp), J)
  for (b in 1:J) P[[b]] <- a[[b]] <- astar[[b]] <- matrix(NA, pjs[[b]], N + 1)
  for (b in 1:J) Y[[b]] <- matrix(NA, nb_ind, N + 1)
  for (n in 1:N) 
  {
     if (verbose) 
      cat(paste0("Computation of the RGCCA block components #", n, " is under progress...\n"))
    if (is.vector(tau))  
      rgcca.result <- rgccak(R, C, tau = tau, scheme = scheme,init = init, bias = bias, tol = tol, verbose = verbose,na.rm=na.rm)
    else rgcca.result <- rgccak(R, C, tau = tau[n, ], scheme = scheme, init = init, bias = bias, tol = tol, verbose = verbose,na.rm=na.rm)
    if (!is.numeric(tau)) 
      tau_mat[n, ] = rgcca.result$call$tau
    AVE_inner[n] <- rgcca.result$AVE_inner
    crit[[n]] <- rgcca.result$crit
    # deflation
    for (b in 1:J) Y[[b]][, n] <- rgcca.result$Y[, b]
    defla.result <- defl.select(rgcca.result$Y, R, ndefl , n, nbloc = J)
    R <- defla.result$resdefl
    for (b in 1:J) P[[b]][, n] <- defla.result$pdefl[[b]]
    for (b in 1:J) a[[b]][, n] <- rgcca.result$a[[b]]
    if (n == 1)
    {
      for (b in 1:J) astar[[b]][, n] <- rgcca.result$a[[b]] 
    }
    else {
      for (b in 1:J)
      {
        astar[[b]][, n] <- rgcca.result$a[[b]] -  astar[[b]][, (1:n - 1), drop = F] %*% drop(t(a[[b]][,n]) %*% P[[b]][, 1:(n - 1), drop = F])                    
      }
    }
  }
  if (verbose) 
    cat(paste0("Computation of the RGCCA block components #", 
               N + 1, " is under progress ... \n"))
  if (is.vector(tau)) 
    rgcca.result <- rgccak(R, C, tau = tau, scheme = scheme, init = init, bias = bias, tol = tol, verbose = verbose)
  else rgcca.result <- rgccak(R, C, tau = tau[N + 1, ], scheme = scheme, init = init, bias = bias, tol = tol, verbose = verbose)
  crit[[N + 1]] <- rgcca.result$crit
  if (!is.numeric(tau)) 
    tau_mat[N + 1, ] = rgcca.result$call$tau
  AVE_inner[max(ncomp)] <- rgcca.result$AVE_inner
  for (b in 1:J) {
    Y[[b]][, N + 1] <- rgcca.result$Y[, b]
    a[[b]][, N + 1] <- rgcca.result$a[[b]]
    # ajout
    astar[[b]][, N + 1] <- rgcca.result$a[[b]] - astar[[b]][, (1:N), drop = F] %*% drop(t(a[[b]][, (N + 1)]) %*%P[[b]][, 1:(N), drop = F])
    #astar[[b]][,N+1] <- rgcca.result$a[[b]] - astar[[b]][,(1:N),drop=F] %*% drop( t(a[[b]][,(N+1)]) %*% P[[b]][,1:(N),drop=F] ) 
  
    rownames(a[[b]]) = rownames(astar[[b]]) = colnames(A[[b]])
    rownames(Y[[b]]) = rownames(A[[b]])
    colnames(Y[[b]]) = paste0("comp", 1:max(ncomp))
  }

  for (j in 1:J)
  {
       # AVE_X[[j]] = apply(cor(A[[j]], Y[[j]],use="pairwise.complete.obs")^2, 	2, mean,na.rm=TRUE)
       # print("AVEinner 1")
       # print( AVE_X[[j]])
       # print("AVEinner 2")
       # AVEinner2=diag(cov(Y[[j]]))/sum(diag(cov(A[[j]] )))
       # print(AVEinner2)    
        AVE_X[[j]]=diag(cov2(Y[[j]]))/sum(diag(cov2(A[[j]] )),na.rm=TRUE)
  }
  outer = matrix(unlist(AVE_X), nrow = max(ncomp))
  
  for (j in 1:max(ncomp)) AVE_outer[j] <- sum(pjs * outer[j,],na.rm=na.rm)/sum(pjs)
  Y = shave.matlist(Y, ncomp)
  names(Y)=names(A)
  names(a)=names(A)
  AVE_X = shave.veclist(AVE_X, ncomp)
  AVE <- list(AVE_X = AVE_X, AVE_outer_model = AVE_outer, AVE_inner_model = AVE_inner)
  out <- list(Y = shave.matlist(Y, ncomp), a = shave.matlist(a,ncomp), astar = shave.matlist(astar, ncomp),  tau = tau_mat,
                crit = crit, primal_dual = primal_dual,	AVE = AVE,A=A,call=call)

   class(out) <- "rgccad"

  return(out)
}
