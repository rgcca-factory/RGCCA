# Product for Matrix (pm) is a generalization of the matricial product %*% for matrices with missing data. Missing data are replaced by 0.
# @param M1  A matrix with n1 lines and p columns
# @param M2  A matrix with p lines and n2 columns
# @param na.rm if TRUE calculates the matricial product only on available data. Else returns NA.
# @return \item{X}{The resulting matrix with n1 lines and n2 columns}
# @title Product for Matrices with missing data (pm)

pm <- function(M1,M2,na.rm=TRUE){
  if(na.rm)
  {
    M1[is.na(M1)]=0
    M2[is.na(M2)]=0
    return(M1%*%M2)
  }
  else
  {
    return(M1%*%M2)
  }
}
