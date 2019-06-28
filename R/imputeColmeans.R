#' Impute a list of matrices by the colmeans
#' @param A A list of matrices
#' @return List of the matrices imputed by colmeans
#' @title Product for Matrices with missing data (pm)
#' @examples
#' #############
#' # Example 1 #
#' #############
#' A=matrix(rnorm(15),3,5);A[3,5]=NA
#' B=matrix(rnorm(15),3,5);A[3,]=NA
#' imputeColmeans(list(A,B))
#' @export imputeColmeans
imputeColmeans=function(A)
{
	imputeAvg=function(x){x[is.na(x)]=mean(x,na.rm=T);return(x)}
	Aimp=lapply(A,function(M){return(apply(M,2,"imputeAvg"))})
	return(Aimp)
}