#' imputeRGCCA allows to choose the imputation method before running RGCCA
#' @inheritParams select_analysis
#' @param blocks A list that contains the J blocks of variables X1, X2, ..., XJ. 
#' Block j is a matrix of dimension $n x p_j$ where $p_j$ is the number of 
#' variables in X_j.
#' @param method  Either a character corresponding to the used method 
#' ("complete","knn","em","sem") or a function taking a list of J blocks (A) as 
#' only parameter and returning the imputed list. 
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
#' @param prescaling A logical value indicating if the scaling should be done outside of the function.
#' @return \item{Y}{A list of \eqn{J} elements. Each element of \eqn{Y} is a matrix that contains the analysis components for the corresponding block.}
#' @return \item{a}{A list of \eqn{J} elements. Each element of \eqn{a} is a matrix that contains the outer weight vectors for each block.}
#' @return \item{astar}{A list of \eqn{J} elements. Each element of astar is a matrix defined as Y[[j]][, h] = A[[j]]\%*\%astar[[j]][, h].}
#' @return \item{C}{A symmetric matrix (J*J) that describes the relationships between blocks}
#' @return \item{tau}{Either a 1*J vector or a \eqn{\mathrm{max}(ncomp) \times J} matrix containing the values 
#' of the regularization parameters (default: tau = 1, for each block and each dimension). Tau varies from 0 (maximizing the correlation) to 1 (maximizing the covariance).
#' If tau = "optimal" the regularization paramaters are estimated for each block and each dimension using the Schafer and Strimmer (2005)
#' analytical formula . If tau is a \eqn{1\times J} vector, tau[j] is identical across the dimensions of block \eqn{\mathbf{X}_j}. 
#' If tau is a matrix, tau[k, j] is associated with \eqn{\mathbf{X}_{jk}} (\eqn{k}th residual matrix for block \eqn{j}). It can be estimated by using \link{rgcca_permutation}.}
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

rgccaNa=function (blocks,method, connection = 1 - diag(length(A)), tau = rep(1, length(A)),    ncomp = rep(1, length(A)), scheme = "centroid", scale = TRUE,   init = "svd", bias = TRUE, tol = 1e-08, verbose = TRUE,
                  scale_block=TRUE,prescaling=FALSE,quiet=FALSE)
{ 
    A=blocks
    C=connection
    call=list(A=A,method=method, C =C, tau = tau,    ncomp = ncomp, scheme = scheme, scale = scale,   init = init, bias = bias, tol =tol, verbose = verbose,scale_block=scale_block)

  nvar = sapply(A, NCOL)

  shave.matlist <- function(mat_list, nb_cols) mapply(function(m,nbcomp) m[, 1:nbcomp, drop = FALSE], mat_list, nb_cols, SIMPLIFY = FALSE)
    shave.veclist <- function(vec_list, nb_elts) mapply(function(m, nbcomp) m[1:nbcomp], vec_list, nb_elts, SIMPLIFY = FALSE)
	A0=A
	indNA=lapply(A,function(x){return(which(is.na(x),arr.ind=TRUE))})
  na.rm=FALSE
  if(method=="complete"){A2=intersection_list(A)}

  if(method=="nipals"){na.rm=TRUE;A2=A}

   resRgcca=rgccad(A2,C=C,ncomp=ncomp,verbose=verbose,scale=scale,init=init,scale_block=scale_block,tau=tau,scheme=scheme,tol=tol,prescaling=prescaling,quiet=quiet)
  out=list(imputedA=A2,rgcca=resRgcca,method,indNA=indNA)
	return(out)

}
