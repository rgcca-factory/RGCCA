#' The function initsvd() is called by rgcca() and does not have to be used by the user. 
#' initsvd() enables the computation of initial scores of subjects for RGCCA based on SVD decomposition
#' If missing values, they are imputed by colmeans
#' @param X  A matrix with n lines and p columns
#' @return A matrix with n lines and n columns
#' @title Initialisation by SVD decomposition of X
#' @export initsvd2
initsvd <- function(X) {
# verifier le scale
  n = NROW(X)
  p = NCOL(X)
  #rowNa=apply(is.na(X),1,sum,na.rm=T)==dim(X)[2]
  vecMoyenne=apply(X,2,mean,na.rm=T)
  matMoyenne=matrix(rep(vecMoyenne,n),n,p)
  X[is.na(X)]=matMoyenne[is.na(X)]
  nbNa=sum(is.na(X))
  if(nbNa>0)
  {
	  res=nipals(X, ncomp = 1, center = FALSE, scale = FALSE)
	  ifelse(n>=p, return(-as.vector(res$loadings)), return(-as.vector(res$scores)) )
	}
  if(nbNa==0){ifelse(n>=p, return(svd(X,nu=0,nv=1)$v), return(svd(X,nu=1,nv=0)$u) )}
}
