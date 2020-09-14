#' imputeRGCCA allows to choose the imputation method before running RGCCA
#' @inheritParams select_analysis
#' @param blocks A list of matrices giving the \eqn{J} blocks of variables \eqn{\mathbf{X_1}, \mathbf{X_2}, ..., \mathbf{X_J}}.
#' @param method  Either a character corresponding to the used method ("complete","knn","em","sem") or a function taking a list of J blocks (A) as only parameter and returning the imputed list. 
#' \itemize{
#' \item{\code{"mean"}}{ corresponds to an imputation by the colmeans}
#' \item{\code{"complete"}}{ corresponds to run RGCCA only on the complete subjects (subjects with missing data are removed)}
#' \item{\code{"nipals"}}{ corresponds to run RGCCA on all available data (NIPALS algorithm)}
#' \item{\code{"em"}}{ corresponds to impute the data with EM-type algorithms}
#' \item{\code{"sem"}}{ corresponds to impute the data with EM-type algorithms with superblock approach}
#' \item{\code{"knn1"}}{ corresponds to impute the data with the 1-Nearest Neighbor. 1 can be replace by another number (such as knn3) to impute with the 3-Nearest Neighbors.}}
#' @param tau Either a 1*J vector or a \eqn{\mathrm{max}(ncomp) \times J} matrix containing the values 
#' of the regularization parameters (default: tau = 1, for each block and each dimension). Tau varies from 0 (maximizing the correlation) to 1 (maximizing the covariance).
#' If tau = "optimal" the regularization paramaters are estimated for each block and each dimension using the Schafer and Strimmer (2005)
#' analytical formula . If tau is a \eqn{1\times J} vector, tau[j] is identical across the dimensions of block \eqn{\mathbf{X}_j}. 
#' If tau is a matrix, tau[k, j] is associated with \eqn{\mathbf{X}_{jk}} (\eqn{k}th residual matrix for block \eqn{j}). It can be estimated by using \link{rgcca_permutation}.
#' @param scale A logical value indicating if each block is normalised and divided by the square root of its number of variables and then divided by the square root of its number of variables.
#' @param scale_block A logical value indicating if each block have the same weight in the RGCCA analysis. Otherwise, the weight of each block depends on the number of variables of the block
#' @param verbose  A logical value indicating if the progress of the analysis will be reported while computing.
#' @param quiet A logical value indicating if it should not print warnings
#' @param init A character giving the mode of initialization to use in the algorithm. The alternatives are either by Singular Value Decompostion ("svd") or random ("random") (default: "svd").
#' @param bias A logical value for biaised (\eqn{1/n}) or unbiaised (\eqn{1/(n-1)}) estimator of the var/cov (default: bias = TRUE).
#' @param tol An integer giving the value for stopping the algorithm convergence.
#' @param knn.k  Used only if missing values in the blocks are estimated by k-NN methods. An integer giving the number of k nearest neighbors. Can also be "auto" for automatic selection.
#' @param knn.output A character among "mean", "random" or "weightedMean" : Used only if missing values in the blocks are estimated by k-NN methods. Returns respectively the average of the k nearest neigbors, one selected randomly, or an average weighted by the distance of the k NN
#' @param knn.klim Used only if missing values in the blocks are estimated by k-NN methods, and if knn.k is "auto". An integer giving the k limits (if k is not a number, optimal k between klim[1] and klim[2] is calculated )
#' @param knn.scale_block Used only if missing values in the blocks are estimated by k-NN methods. A logical value indicating if the distance for Nearest Neigbors takes the size of blocks into account
#' @param pca.ncp An integer giving the number of components chosen in PCA 
#' @param prescaling A logical value indicating if the scaling should be done outside of the function.
#' @param ni An integer giving the number of iterations for em or sem methods
#' @return \item{Y}{A list of \eqn{J} elements. Each element of \eqn{Y} is a matrix that contains the RGCCA components for the corresponding block.}
#' @return \item{a}{A list of \eqn{J} elements. Each element of \eqn{a} is a matrix that contains the outer weight vectors for each block.}
#' @return \item{astar}{A list of \eqn{J} elements. Each element of astar is a matrix defined as Y[[j]][, h] = A[[j]]\%*\%astar[[j]][, h].}
#' @return \item{C}{A symmetric matrix (J*J) that describes the relationships between blocks}
#' @return \item{tau}{A vector or matrix that contains the values of the shrinkage parameters applied to each block and each dimension (user specified).}
#' @return \item{scheme}{A character or a function giving the link function for 
#' covariance maximization among "horst" (the identity function), "factorial"
#'  (the squared values), "centroid" (the absolute values). Only, the horst 
#'  scheme penalizes structural negative correlation. The factorial scheme 
#'  discriminates more strongly the blocks than the centroid one.}
#' @return \item{ncomp}{A vector of 1*J integers giving the number of component for each blocks}
#' @return \item{crit}{A vector that contains the values of the criteria across iterations.}
#' @return \item{mode}{A \eqn{1 \times J} vector that contains the formulation ("primal" or "dual") applied to each of the \eqn{J} blocks within the RGCCA alogrithm} 
#' @return \item{AVE}{indicators of model quality based on the Average Variance Explained (AVE): AVE(for one block), AVE(outer model), AVE(inner model).}
#' @references Tenenhaus A. and Tenenhaus M., (2011), Regularized Generalized Canonical Correlation Analysis, Psychometrika, Vol. 76, Nr 2, pp 257-284.
#' @references Tenenhaus A. et al., (2013), Kernel Generalized Canonical Correlation Analysis, submitted.
#' @references Schafer J. and Strimmer K., (2005), A shrinkage approach to large-scale covariance matrix estimation and implications for functional genomics. Statist. Appl. Genet. Mol. Biol. 4:32.
#' @title Regularized Generalized Canonical Correlation Analysis (RGCCA) 
#' @export
#' @examples
#' data(Russett)
#' X_agric =as.matrix(Russett[,c("gini","farm","rent")])
#' X_ind = as.matrix(Russett[,c("gnpr","labo")])
#' X_polit = as.matrix(Russett[ , c("demostab", "dictator")])
#' X_agric[c(2,4),]=NA
#' X_ind[1,]=NA
#' X_polit[5,1]=NA
#' A = list(agri=X_agric, ind=X_ind, polit=X_polit)
#' rgccaNa(A,method="nipals")
#' rgccaNa(A,method="knn2")

rgccaNa=function (blocks,method, connection = 1 - diag(length(A)), tau = rep(1, length(A)),    ncomp = rep(1, length(A)), scheme = "centroid", scale = TRUE,   init = "svd", bias = TRUE, tol = 1e-08, verbose = TRUE,
                  scale_block=TRUE,knn.k="all",knn.output="weightedMean",knn.klim=NULL,knn.scale_block=TRUE,pca.ncp=1,ni=20,prescaling=FALSE,quiet=FALSE)
{ 
  #  call=match.call() 
    A=blocks
    C=connection
    call=list(A=A,method=method, C =C, tau = tau,    ncomp = ncomp, scheme = scheme, scale = scale,   init = init, bias = bias, tol =tol, verbose = verbose,scale_block=scale_block,knn.k=knn.k,knn.output=knn.output,knn.klim=knn.klim,knn.scale_block=scale_block,pca.ncp=pca.ncp)

  nvar = sapply(A, NCOL)
  superblockAsList=function(superblock,A)
  {
    Alist=list()
    nvar = sapply(A, NCOL)
    for(j in 1:length(nvar))
    {
      if(j==1){sel=1:nvar[1]}else{debut=sum(nvar[1:(j-1)])+1;fin=debut+(nvar[j]-1);sel=debut:fin}
      Alist[[j]]=as.matrix(superblock[,sel])
      colnames( Alist[[j]])=colnames(A[[j]])
    }
    names(Alist)=names(A)
    return(Alist)
  }
  shave.matlist <- function(mat_list, nb_cols) mapply(function(m,nbcomp) m[, 1:nbcomp, drop = FALSE], mat_list, nb_cols, SIMPLIFY = FALSE)
    shave.veclist <- function(vec_list, nb_elts) mapply(function(m, nbcomp) m[1:nbcomp], vec_list, nb_elts, SIMPLIFY = FALSE)
	A0=A
	indNA=lapply(A,function(x){return(which(is.na(x),arr.ind=TRUE))})
  na.rm=FALSE
  if(method=="complete"){A2=intersection_list(A)}
  if(method=="mean"){		 A2=imputeColmeans(A) }
  if(is.function(method))
  {
      A2=method(A)
  }
# 	if(method=="pca")	
# 	{  
# 	  imputedSuperblock= imputePCA(X=do.call(cbind,A), ncp = pca.ncp, scale = TRUE, method ="em")$completeObs 
# 	  A2=superblockAsList(imputedSuperblock, A)
# 	}
# 
#   if(method=="rpca")
#   {
#     imputedSuperblock= missMDA::imputePCA(do.call(cbind,A), ncp = pca.ncp, scale = TRUE, method ="regularized")$completeObs 
#      A2=superblockAsList(imputedSuperblock, A)
#   }   
# #
# #	if(method=="rgccaPca"){	  A2= imputeSuperblock(A,method="em",opt="rgcca",ncp=ncp,scaleBlock=scaleBlock)}
# 	if(method=="mfa")	
# 	{	 
# 	  imputedSuperblock=imputeMFA(X=do.call(cbind,A), group=nvar, ncp = 1, type=rep("s",length(nvar)), method = "em")$completeObs
# 	  A2=superblockAsList(imputedSuperblock, A)
# 	}

 	if(method=="iterativeSB")	{	  A2=imputeSB(A=A,ncomp=ncomp,scale=scale,scale_block=scale_block,tau=tau,tol=tol,ni=ni)$A	}
    if(method=="em")	{	  A2=imputeEM(A=A,ncomp=ncomp,scale=scale,scale_block=scale_block,tau=tau,naxis=1,ni=ni,C=C,tol=tol,verbose=verbose,reg="y",quiet=quiet)$A	}
   if(substr(method,1,3)=="sem")
   {
     if(substr(method,4,4)=="")
     {
       A2=imputeEM(A=A,superblock=TRUE,ncomp=ncomp,scale=scale,scale_block=scale_block,tau=tau,naxis=1,ni=ni,C=C,tol=tol,verbose=FALSE,reg="y",quiet=quiet)$A
     }
     else
     {
       A2=imputeEM(A=A,superblock=TRUE,ncomp=ncomp,scale=scale,scale_block=scale_block,tau=tau,naxis=as.numeric(substr(method,4,4)),ni=ni,C=C,tol=tol,verbose=FALSE,reg="y",quiet=quiet)$A
     }
   }
 # if(method=="old"){}
  if(method=="emo")	{	  A2=imputeEM(A=A,ncomp=ncomp,scale=scale,scale_block=scale_block,tau=tau,naxis=1,ni=ni,C=C,tol=tol,verbose=FALSE,reg="no")$A	}
  if(method=="emw")	{	  A2=imputeEM(A=A,ncomp=ncomp,scale=scale,scale_block=scale_block,tau=tau,naxis=1,ni=ni,C=C,tol=tol,verbose=FALSE,reg="w")$A	}
#  if(method=="semy")	{	  A2=imputeEM(A=A,ncomp=ncomp,superblock=TRUE,scale=scale,scale_block=scale_block,tau=tau,naxis=1,ni=50,C=C,tol=tol,verbose=verbose,reg="y")$A[1:length(A)]	}
#  if(method=="semw")	{	  A2=imputeEM(A=A,ncomp=ncomp,superblock=TRUE,scale=scale,scale_block=scale_block,tau=tau,naxis=1,ni=50,C=C,tol=tol,verbose=verbose,reg="w")$A[1:length(A)]	}
  
  
  if(method=="nipals"){na.rm=TRUE;A2=A}
  
  if(substr(method,1,3)=="knn")
  {
      if(substr(method,4,4)=="A")
      {
        A2=imputeNN(A ,output=knn.output,k="all",klim=knn.klim,scale_block=knn.scale_block);method=paste(method,":",knn.k,sep="")
      }
      else
      {
        A2=imputeNN(A ,output=knn.output,k=as.numeric(substr(method,4,4)),klim=knn.klim,scale_block=knn.scale_block);method=paste(method,":",knn.k,sep="")
      }
  }
  if(method!="imputeInRgcca1"&&method!="imputeInRgcca2"&&method!="imputeInRgccaSB"&&method!="imputeInRgccaLL"){resRgcca=rgccad(A2,C=C,ncomp=ncomp,verbose=verbose,scale=scale,init=init,scale_block=scale_block,tau=tau,scheme=scheme,tol=tol,estimateNA="no",prescaling=prescaling,quiet=quiet)}
  if(method=="imputeInRgcca1"){resRgcca=rgccad(A,C=C,ncomp=ncomp,init=init,verbose=FALSE,scale=scale,scale_block=scale_block,tau=tau,scheme=scheme,tol=tol,estimateNA="iterative",prescaling=prescaling,quiet=quiet);A2=resRgcca$imputedA;}
  if(method=="imputeInRgcca2"){resRgcca=rgccad(A,C=C,ncomp=ncomp,init=init,verbose=FALSE,scale=scale,scale_block=scale_block,tau=tau,scheme=scheme,tol=tol,estimateNA="first",prescaling=prescaling,quiet=quiet);A2=resRgcca$imputedA;}
  if(method=="imputeInRgccaSB"){resRgcca=rgccad(A,C=C,ncomp=ncomp,init=init,verbose=FALSE,scale=scale,scale_block=scale_block,tau=tau,scheme=scheme,tol=tol,estimateNA="superblock",prescaling=prescaling,quiet=quiet);A2=resRgcca$imputedA[1:length(A)];}
  if(method=="imputeInRgccaLL"){resRgcca=rgccad(A,C=C,ncomp=ncomp,init=init,verbose=TRUE,scale=scale,scale_block=scale_block,tau=tau,scheme=scheme,tol=tol,estimateNA="lebrusquet",prescaling=prescaling,quiet=quiet);A2=resRgcca$imputedA[1:length(A)];}
  out=list(imputedA=A2,rgcca=resRgcca,method,indNA=indNA)
	return(out)

}
