#' determine_patternNA
#' Determines the pattern of missing values and allows to use it to create datasets with the same pattern.
#' @param A list of blocks
#' @param graph if TRUE, the pattern of missing values is plotted
#' @return a list containing the percentage of missing values per variable (pctNA), the percentage of missing values per block (pctNAbyBlock), and the complete individuals percentage (completeSubjectByBlock)
#' @examples 
#' X1=matrix(rnorm(150),30,5)
#'X2=matrix(rnorm(150),30,5)
#'X3=matrix(rnorm(150),30,5)
#'X1[1:10,1]=NA
#'X1[2,2]=NA
#'X1[11,]=NA
#'X2[1:10,1]=NA
#'colnames(X1)=paste0("A",1:5)
#'colnames(X2)=paste0("B",1:5)
#'colnames(X3)=paste0("C",1:5)
#'A=list(bloc1=X1,bloc2=X2,bloc3=X3)
#'determine_patternNA(A)
#' @export
#'@importFrom graphics image
#'@importFrom stats na.omit
determine_patternNA=function(A,graph=TRUE)
{
    # 

    #A=checkSize(A)
    #A=checkRownames(A)
    A=check_blocks(A,add_NAlines = TRUE)
    pctNA=lapply(A,function(x)
        {res=apply(x,2,function(y){
            return(sum(is.na(y))/length(y))
        })
        return(res)})
    pctNAbyBlock=sapply(A,function(X){return(sum(is.na(X))/(dim(X)[1]*dim(X)[2]))})
    completeSubjectByBlock=sapply(A,function(x)
        {
            res=apply(x,1,function(y)
                {
                    return(sum(is.na(y))==0)
                })
        })
  if(graph)
  {
      nvar=sapply(A,NCOL)
      par(mfrow=c(1,length(A)))
      par(mar=c(2,1,4,1))
      for(i in 1:length(A))
      {
          mat=apply(is.na(A[[i]]),2,rev)
          image(t(mat),main=paste0(names(A)[i],"\n(",nvar[i]," var.,",sum(completeSubjectByBlock[,i]), "/",NROW(A[[i]])," complete ind.)"),xaxt="n",yaxt="n",col=c("light blue","black"))
      }
  }
  return(list(pctNA=pctNA,pctNAbyBlock=pctNAbyBlock,completeSubjectByBlock=completeSubjectByBlock))
}
