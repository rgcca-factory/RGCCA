determine_patternNA=function(A,graph=TRUE)
{
    # 

    #A=checkSize(A)
    #A=checkRownames(A)
    
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
          image(t(is.na(A[[i]])),main=paste0(names(A)[i],"\n(",nvar[i]," var.,",sum(completeSubjectByBlock[,i]), "/",NROW(A[[i]])," complete ind.)"),xaxt="n",yaxt="n",col=c("light blue","black"))
      }
  }
  return(list(pctNA=pctNA,pctNAbyBlock=pctNAbyBlock,completeSubjectByBlock))
}
