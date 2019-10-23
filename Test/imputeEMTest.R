# '# ImputeEM test
# 
# '''
#setwd("/home/caroline.peltier/Bureau/RGCCA")
# --enter here the name of the git RGCCA working directory 
#library(RGCCA)
library(MASS)
library(nipals)

namesFiles=dir("./R")
sapply(namesFiles,function(x){source(paste0("./R/",x))})

# Russetts
X_agric =as.matrix(Russett[,c("gini","farm","rent")])
X_ind = as.matrix(Russett[,c("gnpr","labo")])
X_polit = as.matrix(Russett[ , c("demostab", "dictator")])
X_agric[c(2,4),]=NA
X_ind[1,]=NA
X_polit[5,1]=NA
A = list(agri=X_agric, ind=X_ind, polit=X_polit)
A_ref=list(agri=as.matrix(Russett[,c("gini","farm","rent")]),ind=as.matrix(Russett[,c("gnpr","labo")]),polit=as.matrix(Russett[ , colnames(Russett)%in%c("demostab", "dictator")]))
A_ref2=lapply(A_ref,scale)
A2=lapply(A,scale)
# pour 1 axe, non superblock
testDataEM=imputeEM(A=A,ncomp=rep(1,3),scale=TRUE,sameBlockWeight=TRUE,tau=rep(1,3),naxis=1,ni=50,C=matrix(1,3,3)-diag(3),tol=1e-8,scheme="centroid",verbose=TRUE)
testDataEM2=imputeEM(A=A,ncomp=rep(1,3),scale=FALSE,sameBlockWeight=TRUE,tau=rep(1,3),naxis=1,ni=50,C=matrix(1,3,3)-diag(3),tol=1e-8,scheme="centroid",verbose=TRUE)
testDataEM3=imputeEM(A=A,ncomp=rep(1,3),scale=TRUE,sameBlockWeight=FALSE,tau=rep(1,3),naxis=1,ni=50,C=matrix(1,3,3)-diag(3),tol=1e-8,scheme="centroid",verbose=TRUE)
testDataEM4=imputeEM(A=A,ncomp=rep(1,3),scale=FALSE,sameBlockWeight=FALSE,tau=rep(1,3),naxis=1,ni=50,C=matrix(1,3,3)-diag(3),tol=1e-8,scheme="centroid",verbose=TRUE)

#testDataEMSB=imputeEM(A=A2,superblock = TRUE,ncomp=rep(1,3),scale=TRUE,sameBlockWeight=TRUE,tau=rep(1,3),naxis=1,ni=50,C=matrix(1,3,3)-diag(3),tol=1e-6,scheme="centroid")
testDataEM=imputeEM(A=A,superblock=TRUE,ncomp=rep(1,3),scale=TRUE,sameBlockWeight=TRUE,tau=rep(1,3),naxis=1,ni=50,C=matrix(1,3,3)-diag(3),tol=1e-8,scheme="centroid",verbose=TRUE)


plot(testDataEM$crit,pch=16,main="RGCCA criterion")
plot(testDataEM$stab,pch=16,main="Stability")
plot(testDataEM$obj,pch=16,main="RMSE")

# pour 2 comp, 1 axe, pas superblock
testDataEM=imputeEM(A=A,noise=TRUE,ncomp=rep(2,3),scale=TRUE,sameBlockWeight=TRUE,tau=rep(1,3),naxis=1,ni=50,C=matrix(1,3,3)-diag(3),tol=1e-6,scheme="centroid")
testDataEM2=imputeEM(A=A,noise=TRUE,ncomp=rep(2,3),scale=TRUE,sameBlockWeight=TRUE,tau=rep(2,3),naxis=1,ni=50,C=matrix(1,3,3)-diag(3),tol=1e-6,scheme="centroid")


# pour 1 axe, superblock
testDataEMSB=imputeEM(A=A,superblock=TRUE,ncomp=rep(1,3),scale=TRUE,sameBlockWeight=TRUE,tau=rep(1,3),naxis=1,ni=50,C=matrix(1,3,3)-diag(3),tol=1e-6,scheme="centroid")
testDataEMSB=imputeEM(A=A2,superblock=TRUE,ncomp=rep(1,3),scale=TRUE,sameBlockWeight=TRUE,tau=rep(1,3),naxis=1,ni=50,C=matrix(1,3,3)-diag(3),tol=1e-6,scheme="centroid")


plot(testDataEMSB$crit,pch=16,main="RGCCA criterion")
plot(testDataEMSB$stab,pch=16,main="Stability",ylim=c(0,1))
plot(testDataEMSB$obj,pch=16,main="RMSE")

# pour superblock plusieurs axes
testDataEMSB=imputeEM(A=A,noise=TRUE,superblock=TRUE,ncomp=rep(2,3),scale=TRUE,sameBlockWeight=FALSE,tau=rep(1,3),naxis=1,ni=50,C=matrix(1,3,3)-diag(3),tol=1e-6,scheme="centroid")
testDataEMSB2=imputeEM(A=A,ncomp=rep(2,3),superblock=TRUE,scale=TRUE,sameBlockWeight=TRUE,tau=rep(1,3),naxis=2,ni=50,C=matrix(1,3,3)-diag(3),tol=1e-6,scheme="centroid")
testDataEMSB3=imputeEM(A=A,ncomp=rep(1,3),superblock=TRUE,scale=TRUE,sameBlockWeight=TRUE,tau=rep(1,3),naxis=3,ni=50,C=matrix(1,3,3)-diag(3),tol=1e-6,scheme="centroid")



testDataEM$A$ind[1,]
testDataEMSB$A$ind[1,]
testDataEMSB2$A$ind[1,]
testDataEMSB3$A$ind[1,]
as.matrix(Russett[,c("gnpr","labo")])[1,]

testDataEM$A$polit[5,1]
testDataEMSB$A$polit[5,1]
testDataEMSB2$A$polit[5,1]
testDataEMSB3$A$polit[5,1]
as.matrix(Russett[ , c("demostab", "dictator")])[5,1]

testDataEM$A$agri[c(2,4),]
testDataEMSB$A$agri[c(2,4),]
testDataEMSB2$A$agri[c(2,4),]
testDataEMSB3$A$agri[c(2,4),]
as.matrix(Russett[,c("gini","farm","rent")])[c(2,4),]


# Etude du superbloc dans le cas Biosca


setwd("/home/caroline.peltier/Bureau/EtudeNA/Datasets/Biosca/Reference")
refData=readDataset(c("CLI","MRS","VOL"))
A=refData
A[[1]][1:2,]=NA
A[[1]][3:8,]=NA
A[[1]][9:12,]=NA
lapply(A,dim)
refData[[1]][1:4,]
apply(refData[[1]],2,mean)
A2=lapply(A,function(x){apply(x,2,scale)})
testDataEMSB3=imputeEM(A=A2,ncomp=rep(1,3),superblock=TRUE,scale=TRUE,sameBlockWeight=TRUE,tau=rep(1,3),naxis=1,ni=50,C=matrix(1,3,3)-diag(3),tol=1e-6,scheme="centroid")

listNAdataset=createNA(A=referenceDataset,option="block",pNA=patternNA,nAllRespondants=10,output="list")
# TO DO test avec deux composantes et un axe
#testDataEM=imputeEM(A=A,ncomp=rep(2,3),scale=TRUE,sameBlockWeight=TRUE,tau=rep(1,3),naxis=1,ni=50,C=matrix(1,3,3)-diag(3),tol=1e-8)
#testDataEM=imputeEM(A=A,ncomp=rep(2,3),superblock=TRUE,scale=TRUE,sameBlockWeight=TRUE,tau=rep(1,3),naxis=1,ni=50,C=matrix(1,3,3)-diag(3),tol=1e-8)

# TO DO test avec deux composantes et deux axes
testDataEMSB2=imputeEM(A=A,ncomp=rep(2,3),superblock=TRUE,scale=TRUE,sameBlockWeight=TRUE,tau=rep(1,3),naxis=2,ni=50,C=matrix(1,3,3)-diag(3),tol=1e-6,scheme="centroid")



w=c(5,3,2,1,4,5,6,2,1)
gamma=c(1,-1,3,2,4,5,1,5)
(as.matrix(gamma))%*%t(as.matrix(w))+ rnorm(0,1)
