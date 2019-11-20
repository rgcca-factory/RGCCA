#' The function writeList write 
#' @param A A list of blocks
#' @param wd Adress of the working directory where the files will be created
#' @title write the different blocks of a list in a working directory
#' @examples
#' set.seed(42);X1=matrix(rnorm(35),7,5);X2=matrix(rnorm(28),7,4)
#' A=list(X1,X2)
#' writeList(A,wd=getwd())
#' @export writeList
writeList=function(A,wd=getwd())
{
  setwd(wd)
  for(i in 1:length(A))
  {
    write.table(A[[i]],paste(names(A)[i],".csv",sep=""),sep=";",row.names=TRUE)
  }
}