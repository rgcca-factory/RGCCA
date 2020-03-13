#' imputeRGCCA allows to choose the imputation method before running RGCCA
#' @param blocks  A list that contains the \eqn{J} blocks of variables \eqn{\mathbf{X_1}, \mathbf{X_2}, ..., \mathbf{X_J}}.
#' @param method  Either a character corresponding to the used method ("complete","knn","em","sem") or a function taking a list of J blocks (A) as only parameter and returning the imputed list. 
#' @param connection  A design matrix that describes the relationships between blocks (default: complete design).
#' @param sparsity Used for type="sgcca" or "spls" only. A vector containing the sparsity coefficients (length J, between 0 and 1). It can be estimated by using \link{rgcca_permutation}.
#' @param scheme The value is "horst", "factorial", "centroid" or the g function (default: "centroid").
#' @param scale  If scale = TRUE, each block is standardized to zero means and unit variances (default: TRUE).
#' @param sameBlockWeight TRUE by default : each block have the same weight in the RGCCA analysis. If FALSE, the weight of each block depends on the number of variables of the block
#' @param ncomp  A \eqn{1 \times J} vector that contains the numbers of components for each block (default: rep(1, length(A)), which gives one component per block.). It can be estimated by using \link{rgcca_permutation}.
#' @param verbose  If verbose = TRUE, the progress will be report while computing (default: TRUE).
#' @param quiet If TRUE, does not print warnings
#' @param init The mode of initialization to use in RGCCA algorithm. The alternatives are either by Singular Value Decompostion ("svd") or random ("random") (Default: "svd").
#' @param bias A logical value for biaised or unbiaised estimator of the var/cov (default: bias = TRUE).
#' @param tol The stopping value for convergence.
#' @param knn.k  Used only if missing values in the blocks are estimated by k-NN methods. Number of k nearest neighbors. Can also be "auto" for automatic selection.
#' @param knn.output "mean", "random" or "weightedMean" : Used only if missing values in the blocks are estimated by k-NN methods. Returns respectively the average of the k nearest neigbors, one selected randomly, or an average weighted by the distance of the k NN
#' @param knn.klim Used only if missing values in the blocks are estimated by k-NN methods, and if knn.k is "auto". k limits (if k is not a number, optimal k between klim[1] and klim[2] is calculated )
#' @param knn.sameBlockWeight Used only if missing values in the blocks are estimated by k-NN methods.if TRUE the distance for Nearest Neigbors takes the size of blocks into account
#' @param pca.ncp Number of components chosen in PCA 
#' @param prescaling If TRUE, sgcca does NOT run scaling steps (they were calculated before)
#' @return \item{Y}{A list of \eqn{J} elements. Each element of \eqn{Y} is a matrix that contains the RGCCA components for the corresponding block.}
#' @return \item{a}{A list of \eqn{J} elements. Each element of \eqn{a} is a matrix that contains the outer weight vectors for each block.}
#' @return \item{astar}{A list of \eqn{J} elements. Each element of astar is a matrix defined as Y[[j]][, h] = A[[j]]\%*\%astar[[j]][, h].}
#' @return \item{C}{A design matrix that describes the relation between blocks (user specified).}
#' @return \item{scheme}{The scheme chosen by the user (user specified).}
#' @return \item{ncomp}{A \eqn{1 \times J} vector that contains the numbers of components for each block (user specified).}
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

sgccaNa=function (blocks,method, connection = 1 - diag(length(A)), sparsity = rep(1, length(A)),    ncomp = rep(1, length(A)), scheme = "centroid", scale = TRUE,   init = "svd", bias = TRUE, tol = 1e-08, verbose = TRUE,sameBlockWeight=TRUE,knn.k="all",knn.output="weightedMean",knn.klim=NULL,knn.sameBlockWeight=TRUE,pca.ncp=1,prescaling=FALSE,quiet=FALSE)
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
  if(method=="complete"){A2=intersection(A)}
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
        A2=imputeNN(A ,output=knn.output,k="all",klim=knn.klim,sameBlockWeight=knn.sameBlockWeight);method=paste(method,":",knn.k,sep="")
      }
      else
      {
        A2=imputeNN(A ,output=knn.output,k=as.numeric(substr(method,4,4)),klim=knn.klim,sameBlockWeight=knn.sameBlockWeight);method=paste(method,":",knn.k,sep="")
      }
  }

  resRgcca=sgcca(A2,sparsity=sparsity,ncomp=ncomp,verbose=FALSE,scale=scale,sameBlockWeight=sameBlockWeight,scheme=scheme,tol=tol,prescaling=prescaling,quiet=quiet)
 return(list(imputedA=A2,rgcca=resRgcca,method,indNA=indNA))

}
