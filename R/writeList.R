#' The function writeList write 
#' @param A A list of blocks
#' @param wd Adress of the working directory where the files will be created
#' @return 
#' @title write the different blocks of a list in a working directory
#' @examples
#' @export writeList
writeList=function(A,wd=getwd())
{
  setwd(wd)
  for(i in 1:length(A))
  {
    write.table(A[[i]],paste(names(A)[i],".csv",sep=""),sep=";",row.names=TRUE)
  }
}