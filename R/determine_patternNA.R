#' get_patternNA
#' Determines the pattern of missing values and allows to use it to create datasets with the same pattern.
#' @param blocks list of blocks
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
#'get_patternNA(A)
#' @export
#'@importFrom graphics image
#'@importFrom stats na.omit
#'@seealso \code{\link[RGCCA]{plot.patternNA}}
get_patternNA=function(blocks)
{
    # 

    #A=checkSize(A)
    #A=checkRownames(A)
  #  check_blocks(A,add_NAlines=TRUE)
     blocks=check_blocks(blocks,add_NAlines = TRUE)
    pctNA=lapply(blocks,function(x)
        {
            if(!is.null(dim(x))||dim(x)[2]==1)
            {
                res=apply(x,2,function(y){
                    return(sum(is.na(y))/length(y))
                })
            }
            else
            {
                res=sum(is.na(x))/NROW(x)
                return(res) 
            }
        return(res)
        })
    pctNAbyBlock=sapply(blocks,function(X){return(sum(is.na(X))/(dim(X)[1]*dim(X)[2]))})
    completeSubjectByBlock=sapply(blocks,function(x)
        {
            
            res=apply(x,1,function(y)
                {
                    return(sum(is.na(y))==0)
                })
            
        })
    rownames(completeSubjectByBlock)=rownames(blocks[[1]])
    completeSubjectsBool=apply(Reduce(cbind,blocks),1,function(x){sum(is.na(x))==0})
    completeSubjects=rownames(blocks[[1]])[completeSubjectsBool]
    numberOfMissingBlocksPerSubject=    sapply(rownames(blocks[[1]]),function(k)
    {
        
        y=sum(sapply(blocks,function(x)
        { 
            return(sum(!is.na(x[k,]))==0)
        }))
        ;return(y)
    })

   patternNA=list(pctNA=pctNA,pctNAbyBlock=pctNAbyBlock,completeSubjectByBlock=completeSubjectByBlock,completeSubjects=completeSubjects, numberOfMissingBlocksPerSubject= numberOfMissingBlocksPerSubject,blocks=blocks)
   class(patternNA)="patternNA"
   return(patternNA)
}
