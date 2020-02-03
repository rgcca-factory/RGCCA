#' determine_patternNA
#' Determines the pattern of missing values and allows to use it to create datasets with the same pattern.
#' @param A list of blocks
#' @param graph if TRUE, the pattern of missing values is plotted
#' @param legend if TRUE the legend is plotted
#' @param outlierVisible if FALSE, the outliers will be -2 if negative, 2 if positive
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
determine_patternNA=function(A,graph="all",legend=FALSE,outlierVisible=TRUE,scale=TRUE)
{
    # 

    #A=checkSize(A)
    #A=checkRownames(A)
    A=check_blocks(A,add_NAlines = TRUE)
    pctNA=lapply(A,function(x)
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
    pctNAbyBlock=sapply(A,function(X){return(sum(is.na(X))/(dim(X)[1]*dim(X)[2]))})
    completeSubjectByBlock=sapply(A,function(x)
        {
            
            res=apply(x,1,function(y)
                {
                    return(sum(is.na(y))==0)
                })
            
        })
    names(completeSubjectByBlock)=rownames(A[[1]])
    completeSubjectsBool=apply(Reduce(cbind,A),1,function(x){sum(is.na(x))==0})
    completeSubjects=rownames(A[[1]])[completeSubjectsBool]
    numberOfMissingBlocksPerSubject=    sapply(rownames(A[[1]]),function(k)
    {
        
        y=sum(sapply(A,function(x)
        { 
            return(sum(!is.na(x[k,]))==0)
        }))
        ;return(y)
    })
    if(graph=="all")
  {
      nvar=sapply(A,NCOL)
      par(mfrow=c(1,ifelse(legend,length(A)+1,length(A))))
      par(mar=c(2,1,4,1))
      for(i in 1:length(A))
      {
          print(i)
         
          mat=apply(A[[i]],2,rev)
          mat2=apply(mat,2,function(x) scale(x,scale=scale))
          
           minimum=-2*max(apply(mat2,2,function(x) sd(x,na.rm=TRUE)),na.rm=TRUE)
          maximum=2*max(apply(mat2,2,function(x) sd(x,na.rm=TRUE)),na.rm=TRUE)
          if(!outlierVisible)
          {
              mat2[mat2< minimum]=minimum
              mat2[mat2>maximum]=maximum
          }
         print(minimum)
         print(maximum)
          plot(NULL,xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n",bty="n",main=paste0(names(A)[i],"\n",nvar[i]," var.\n",sum(completeSubjectByBlock[,i]), "/",NROW(A[[i]])," complete ind."))
          rect(0,0,1,1,col="black")
          image(t(mat2),breaks=seq(minimum,maximum,length.out=51),xaxt="n",yaxt="n",col=rainbow(50,start=0,end=0.5),add=TRUE)
      }
  }
  if(graph %in% names(A))
  {
      nvar=NCOL(A[[graph]])
      par(mfrow=c(1,1))
      mat=apply(A[[graph]],2,rev)
      par(bg="black")
      image(bg="black",t(mat),main=paste0(names(A)[graph],"\n(",nvar," var.,",sum(completeSubjectByBlock[,graph]), "/",NROW(A[[graph]])," complete ind.)"),xaxt="n",yaxt="n",col=c("light blue","black"))
      par(bg="white")
      }
 if(legend)
 {
     plot(NULL,xlim=c(0,1),ylim=c(0,50),bty="n",xaxt="n",yaxt="n")
     lapply(1:50,function(i){points(0.5,i,col=rainbow(50,start=0,end=0.5)[i],pch=15)})
     text(0.6, 1, "min");   text(0.6, 50, "max")
     
 }
  return(list(pctNA=pctNA,pctNAbyBlock=pctNAbyBlock,completeSubjectByBlock=completeSubjectByBlock,completeSubjects=completeSubjects, numberOfMissingBlocksPerSubject= numberOfMissingBlocksPerSubject,A=A))
}
