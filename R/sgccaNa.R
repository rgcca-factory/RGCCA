#' imputeRGCCA allows to choose the imputation method before running RGCCA
#' @inheritParams select_analysis
#' @inheritParams rgccaNa
#' @inheritParams sgcca
#' @return \item{Y}{A list of \eqn{J} elements. Each element of \eqn{Y} is a matrix that contains the analysis components for the corresponding block.}
#' @return \item{a}{A list of \eqn{J} elements. Each element of \eqn{a} is a matrix that contains the outer weight vectors for each block.}
#' @return \item{astar}{A list of \eqn{J} elements. Each element of astar is a matrix defined as Y[[j]][, h] = A[[j]]\%*\%astar[[j]][, h].}
#' @return \item{C}{A symmetric matrix (J*J) that describes the relationships between blocks}
#' @return \item{scheme}{A character or a function giving the link function for 
#' covariance maximization among "horst" (the identity function), "factorial"
#'  (the squared values), "centroid" (the absolute values). Only, the horst 
#'  scheme penalizes structural negative correlation. The factorial scheme 
#'  discriminates more strongly the blocks than the centroid one}
#' @return \item{ncomp}{A vector of 1*J integers giving the number of component for each blocks}
#' @return \item{crit}{A vector of integer that contains for each component the values of the analysis criteria across iterations.}
#' @return \item{mode}{A \eqn{1 \times J} vector that contains the formulation ("primal" or "dual") applied to each of the \eqn{J} blocks within the RGCCA alogrithm} 
#' @return \item{AVE}{A list of numerical values giving the indicators of model quality based on the Average Variance Explained (AVE): AVE(for each block), AVE(outer model), AVE(inner model).}
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

sgccaNa=function (blocks,method, connection = 1 - diag(length(A)), sparsity = rep(1, length(A)),    ncomp = rep(1, length(A)), scheme = "centroid", scale = TRUE,   init = "svd", bias = TRUE, tol = 1e-08, verbose = TRUE,scale_block=TRUE,knn.k="all",knn.output="weightedMean",knn.klim=NULL,knn.scale_block=TRUE,pca.ncp=1,prescaling=FALSE,quiet=FALSE)
{ 
  A=blocks
  C=connection
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

 resRgcca=sgcca(A2,C=connection,sparsity=sparsity,ncomp=ncomp,verbose=verbose,scale=scale,scale_block=scale_block,scheme=scheme,tol=tol,prescaling=prescaling,quiet=quiet)
 return(list(imputedA=A2,rgcca=resRgcca,method,indNA=indNA))

}
