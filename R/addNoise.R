addNoise=function(resImputeEM,superblock=FALSE)
{
  sigma=resImputeEM$sigma
  stdev=resImputeEM$stdev
  A=resImputeEM$A
  indNA=resImputeEM$indNA
  moy=resImputeEM$moy
  nvar=sapply(A,NCOL)
  J=length(A)
  concatenedBlocks = Reduce(cbind,Alist)
  eps=list()
  if(length(moy)>1)
  {
    for(j in 1:J)
    {
      eps[[j]]=matrix(rnorm(dim(moy[[j]])[1]*dim(moy[[j]])[2],m=0,sd=sigma[[j]]),nrow=dim(moy[[j]])[1],ncol=dim(moy[[j]])[2])
      #superblock=do.call(cbind,A)
      
      #Xhat =(gamma%*%t(w))*stdev + moy
      #Xhat =(gamma%*%t(w)+eps)*stdev + moy=gamma%*%t(w)*stdev+ moy+eps*stdev 
      A[[j]][indNA[[j]]]=A[[j]][indNA[[j]]]+(eps[[j]]*stdev[[j]])[indNA[[j]]]
    }
    
    Alist=list()
    for(j in 1:length(nvar))
    {
      if(j==1){sel=1:nvar[1]}else{debut=sum(nvar[1:(j-1)])+1;fin=debut+(nvar[j]-1);sel=debut:fin}
      Alist[[j]]=as.matrix(concatenedBlocks[,sel])
      colnames( Alist[[j]])=colnames(A[[j]])
    }
    names(Alist)=names(A)   
  }
   
}
