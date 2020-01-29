
 set.seed(42);X1=matrix(rnorm(500),100,5);
 set.seed(22);X2=matrix(rnorm(400),100,4);
 set.seed(2);X3=matrix(rnorm(700),100,7);
 rownames(X1)=rownames(X2)=rownames(X3)=paste("S",1:100,sep="")
 colnames(X1)=paste("A",1:5)
 colnames(X2)=paste("B",1:4)
 colnames(X3)=paste("C",1:7)
 X1[1,]=NA
 X2[7,1]=NA
  X2[5,1]=NA
  X3[3,1:2]=NA
 A=list(X1,X2,X3)
res=MIRGCCA(A,k=3,ni=5,scale=TRUE,sameBlockWeight=TRUE,tau=rep(0,3))
