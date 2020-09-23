# The function initsvd() is called by rgccad() and does not have to be used by the user.
# initsvd() enables the computation of initial scores of subjects for RGCCA based on SVD decomposition
# If missing values, they are imputed by colmeans
# @param X  A matrix with n lines and p columns
# @param dual TRUE by default, allow to study the transposed matrix X when the number of rows is lower that the number of columns
# @return A matrix with n lines and n columns
# @title Initialisation by SVD decomposition of X

initsvd <- function(X,dual=TRUE) {
    # verifier le scale
    n = NROW(X)
    p = NCOL(X)
    #rowNa=apply(is.na(X),1,sum,na.rm=T)==dim(X)[2]
    vecMoyenne=apply(X,2,mean,na.rm=TRUE)
    matMoyenne=matrix(rep(vecMoyenne,n),n,p)
    X[is.na(X)]=matMoyenne[is.na(X)]
    if(dual)
    { 
        ifelse(n>=p,
               {
               return(svd(X,nu=0,nv=1)$v)},
               {
               return(svd(X,nu=1,nv=0)$u)} )
  
    }
    else
    { 
        v=svd(X,nu=0,nv=1)$v
       
        return(v)
   }
}
