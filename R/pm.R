#' Product for Matrix (pm) is a generalization of the matricial product %*% for matrices with missing data. Missing data are replaced by 0.
#' @param M1  A matrix with n1 lines and p columns
#' @param M2  A matrix with p lines and n2 columns
#' @param na.rm if TRUE calculates the matricial product only on available data. Else returns NA.
#' @return \item{X}{The resulting matrix with n1 lines and n2 columns}
#' @title Product for Matrices with missing data (pm)

pm=function(M1,M2,na.rm=TRUE)
{
  M1b=as.matrix(M1); M2b=as.matrix(M2)
  if(dim(M1b)[2]!=dim(M2b)[1]){stop_rgcca("matrices should have a number of rows/columns compatible with matricial product")}
  
  if(na.rm)
  {
    ind1=which(is.na(M1),arr.ind=TRUE)
    ind2=which(is.na(M2),arr.ind=TRUE)
    M1b[ind1]=0
    M2b[ind2]=0
    return(M1b%*%M2b)
  }
  else
  {
    return(M1b%*%M2b)
  }
	
}
