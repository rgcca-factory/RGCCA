addNoise=function(resImputeEM,superblock=FALSE)
{
  sigma=resImputeEM$sigma
  stdev=resImputeEM$stdev
  A=resImputeEM$A
  indNA=resImputeEM$indNA
  moy=resImputeEM$moy
  nvar=sapply(A,NCOL)
  J=length(A)
  concatenedBlocks = Reduce(cbind,A)
  eps=list()
  if(length(moy)>1)
  {
    for(j in 1:J)
    {
      eps[[j]]=matrix(rnorm(dim(A[[j]])[1]*dim(A[[j]])[2],mean=0,sd=sigma[[j]]),nrow=dim(A[[j]])[1],ncol=dim(A[[j]])[2])
      #superblock=do.call(cbind,A)

      #Xhat =(gamma%*%t(w))*stdev + moy
      #Xhat =(gamma%*%t(w)+eps)*stdev + moy=gamma%*%t(w)*stdev+ moy+eps*stdev
      A[[j]][indNA[[j]]]=A[[j]][indNA[[j]]]+(eps[[j]]*stdev[[j]])[indNA[[j]]]
    }
  }
   return(A)
}
