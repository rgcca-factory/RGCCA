#' imputeRGCCA allows to choose the imputation method before running RGCCA
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
#' @return \item{Y}{A list of \eqn{J} elements. Each element of \eqn{Y} is a matrix that contains the RGCCA components for the corresponding block.}
#' @return \item{a}{A list of \eqn{J} elements. Each element of \eqn{a} is a matrix that contains the outer weight vectors for each block.}
#' @return \item{astar}{A list of \eqn{J} elements. Each element of astar is a matrix defined as Y[[j]][, h] = A[[j]]\%*\%astar[[j]][, h].}
#' @return \item{C}{A design matrix that describes the relation between blocks (user specified).}
#' @return \item{tau}{A vector or matrix that contains the values of the shrinkage parameters applied to each block and each dimension (user specified).}
#' @return \item{scheme}{The scheme chosen by the user (user specified).}
#' @return \item{ncomp}{A \eqn{1 \times J} vector that contains the numbers of components for each block (user specified).}
#' @return \item{crit}{A vector that contains the values of the criteria across iterations.}
#' @return \item{mode}{A \eqn{1 \times J} vector that contains the formulation ("primal" or "dual") applied to each of the \eqn{J} blocks within the RGCCA alogrithm} 
#' @return \item{AVE}{indicators of model quality based on the Average Variance Explained (AVE): AVE(for one block), AVE(outer model), AVE(inner model).}
#' @references Tenenhaus A. and Tenenhaus M., (2011), Regularized Generalized Canonical Correlation Analysis, Psychometrika, Vol. 76, Nr 2, pp 257-284.
#' @references Tenenhaus A. et al., (2013), Kernel Generalized Canonical Correlation Analysis, submitted.
#' @references Schafer J. and Strimmer K., (2005), A shrinkage approach to large-scale covariance matrix estimation and implications for functional genomics. Statist. Appl. Genet. Mol. Biol. 4:32.
#' @title Regularized Generalized Canonical Correlation Analysis (RGCCA) 
#' @examples imputeRGCCA(A,C)
#' @export imputeRGCCA

imputeRGCCA=function (A, C = 1 - diag(length(A)), tau = rep(1, length(A)), refData=NULL,    ncomp = rep(1, length(A)), scheme = "centroid", scale = TRUE,   init = "svd", bias = TRUE, tol = 1e-08, verbose = TRUE,na.impute="none",na.niter=10,na.keep=NULL,nboot=10,sameBlockWeight=TRUE,ncp=1,scaleBlock="rgcca") 
{ 

  shave.matlist <- function(mat_list, nb_cols) mapply(function(m,nbcomp) m[, 1:nbcomp, drop = FALSE], mat_list, nb_cols, SIMPLIFY = FALSE)
	shave.veclist <- function(vec_list, nb_elts) mapply(function(m, nbcomp) m[1:nbcomp], vec_list, nb_elts, SIMPLIFY = FALSE)
	A0=A

 if(na.impute=="mean"){		 A2=imputeColmeans(A) }
 if(na.impute=="knn")	{ A2= imputeNN(A)	  }

	if(na.impute=="pca")	{   A2= imputeSuperblock(A,ncp=ncp,method="em",opt="pca") }
	if(na.impute=="rgccaPca"){	  A2= imputeSuperblock(A,method="em",opt="rgcca",ncp=ncp,scaleBlock=scaleBlock)}
	if(na.impute=="mfa")	{	  A2= imputeSuperblock(A,opt="mfa",method="em",scaleBlock=scaleBlock)}
 	if(na.impute=="iterativeSB")	{	  A2=imputeSB(A,tau,niter=na.niter,graph=TRUE,tol=tol,naxis=1)$A	}
	  
	return(A2)

}