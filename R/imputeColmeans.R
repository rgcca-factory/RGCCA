#' Impute a list of matrices by the colmeans
#' @param A A list of matrices
#' @return List of the matrices imputed by colmeans
#' @title Product for Matrices with missing data (pm)

imputeColmeans=function(A)
{
	imputeAvg=function(x){x[is.na(x)]=mean(x,na.rm=T);return(x)}
	Aimp=lapply(A,function(M){return(apply(M,2,"imputeAvg"))})
	return(Aimp)
}