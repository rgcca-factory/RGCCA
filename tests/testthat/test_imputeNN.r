# test with a missing line
 bloc1=matrix(rbind(c(1,1,2),c(1,1.1,2),c(3,5,1),c(2,2,5),c(2,2,5),c(2,2,5),c(1,1,1),c(1,2,2)),8,3);rownames(bloc1)=paste("S",c('21','1','4','6','4','5','a','b'));colnames(bloc1)=c('A','B','C')
 bloc2=matrix(rbind(c(0,3,0),c(1,12,0),c(NA,NA,NA),c(1,1,2),c(1,1,2),c(1,1,2),c(1,1,1),c(1,2,2)),8,3);rownames(bloc2)=paste("S",c('21','1','4','6','4','5','a','b'));colnames(bloc2)=c('D','E','F')
 bloc3=matrix(rbind(c(0,3,1),c(1,NA,0),c(1,1,1),c(1,2,2),c(1,1,1),c(1,2,2),c(1,1,1),c(1,2,2)),8,3);rownames(bloc3)=paste("S",c('21','1','4','6','4','5','a','b'));colnames(bloc3)=c('G','I','H')
 listRes=list(bloc1,bloc2,bloc3)
 res=imputeNN(listRes,k = 1,
              output = "mean",
              klim = NULL,
              scale = TRUE,
              scale_block = TRUE)

 set.seed(42);X1=matrix(rnorm(500),100,5);
  set.seed(22);X2=matrix(rnorm(400),100,4);
  set.seed(2);X3=matrix(rnorm(700),100,7);
  X1[1,]=NA;X2[7,1]=NA;X2[5,1]=NA;X3[3,]=NA;X3[4,]=NA
  rownames(X1)=rownames(X2)=rownames(X3)=paste("S",1:100,sep="")
  colnames(X1)=paste("A",1:5,sep="");colnames(X2)=paste("A",6:9,sep="");
  colnames(X3)=paste("A",10:16,sep="")
 A=list(X1,X2,X3)  
 Ares=imputeNN(A,k=1,output="mean", klim=NULL,scale=TRUE,scale_block=TRUE)
 

 Ares=imputeNN(A,k=1,output="mean", klim=NULL,scale=TRUE,scale_block=TRUE,superblock=TRUE)
 A=list(X2,X2)
 Ares=imputeNN(A,k=1,output="mean", klim=NULL,scale=TRUE,scale_block=TRUE,superblock=FALSE)
 