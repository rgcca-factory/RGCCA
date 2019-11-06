# '# Test MIRGCCA
# 
# '''

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

# test for knn
res=MIRGCCA(A,k=4,ni=5,scale=TRUE,sameBlockWeight=FALSE,tau=rep(0,3),klim=NULL,output="weightedMean",scheme="centroid",tol=1e-16)
res$rgccaList[[2]]$Y[[1]]
res$rgccaList[[3]]$Y[[1]]
plotMIRGCCA(res,multiple="ell",indnames=FALSE)


# test for em without superblock
res=MIRGCCA(A,option="em",C=matrix(1,3,3)-diag(3),k=4,ni=5,scale=TRUE,sameBlockWeight=FALSE,tau=rep(0,3),klim=NULL,output="weightedMean",scheme="centroid",tol=1e-16,superblock=FALSE)
plotMIRGCCA(res,multiple="seg",indnames=FALSE)


# tests sur Russetts
X_agric =as.matrix(Russett[,c("gini","farm","rent")]);
X_ind = as.matrix(Russett[,c("gnpr","labo")]);
X_polit = as.matrix(Russett[ , c("demostab", "dictator")]);
C=matrix(1,3,3)-diag(3)
A = list(X_agric, X_ind, X_polit);
A_miss0=createNA(A=A,option="ponc",pNA=0.03,nAllRespondants=10,output="list")$dat

res=MIRGCCA(A_miss0,option="em",C=matrix(1,3,3)-diag(3),k=4,ni=20,scale=TRUE,sameBlockWeight=FALSE,tau=rep(0,3),klim=NULL,output="weightedMean",scheme="centroid",tol=1e-16,superblock=FALSE)
plotMIRGCCA(res,multiple="seg",indnames=FALSE,cex=0.5)

plotMIRGCCA(res,multiple="ell",indnames=TRUE,cex=0.5)
