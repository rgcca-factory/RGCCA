# '# Test defl.select
# 
# '''
set.seed(42);X1=matrix(rnorm(15),3,5);
set.seed(22);X2=matrix(rnorm(12),3,4);
set.seed(2);X3=matrix(rnorm(12),3,7);
A=list(X1,X2,X3)
yy=cbind(c(1,0,0),c(0,0,1),c(1/sqrt(2),1/sqrt(2),0)) # projection vectors are identity
res=defl.select (yy=yy, rr=A, nncomp=c(1,1,1), nn=1, nbloc=3) 
  
res$pdefl[[1]]==X1[1,]
res$pdefl[[2]]==(X2[3,])

