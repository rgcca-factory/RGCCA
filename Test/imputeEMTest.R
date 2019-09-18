# '# ImputeSB test
# 
# '''
#setwd("/home/caroline.peltier/Bureau/RGCCA")
library(RGCCA)
library(MASS)
library(nipals)
data(Russett)

namesFiles=dir("./R")
sapply(namesFiles,function(x){source(paste0("./R/",x))})



X_agric =as.matrix(Russett[,c("gini","farm","rent")])
X_ind = as.matrix(Russett[,c("gnpr","labo")])
X_polit = as.matrix(Russett[ , c("demostab", "dictator")])
X_agric[c(2,4),]=NA
X_ind[1,]=NA
X_polit[5,1]=NA
A = list(agri=X_agric, ind=X_ind, polit=X_polit)

testDataEM=imputeEM(A=A,ncomp=rep(1,3),scale=TRUE,sameBlockWeight=TRUE,tau=rep(1,3),naxis=1,ni=50,C=matrix(1,3,3)-diag(3),tol=1e-6,scheme="centroid")
plot(testDataEM$crit,pch=16,main="RGCCA criterion")
plot(testDataEM$stab,pch=16,main="Stability")
plot(testDataEM$obj,pch=16,main="RMSE")
testDataEM=imputeEM(A=A,noise=TRUE,ncomp=rep(1,3),scale=TRUE,sameBlockWeight=TRUE,tau=rep(1,3),naxis=1,ni=50,C=matrix(1,3,3)-diag(3),tol=1e-6,scheme="centroid")


testDataEMSB=imputeEM(A=A,superblock=TRUE,ncomp=rep(1,3),scale=TRUE,sameBlockWeight=FALSE,tau=rep(1,3),naxis=1,ni=50,C=matrix(1,3,3)-diag(3),tol=1e-6,scheme="centroid")
testDataEMSB=imputeEM(A=A,noise=TRUE,superblock=TRUE,ncomp=rep(1,3),scale=TRUE,sameBlockWeight=FALSE,tau=rep(1,3),naxis=1,ni=50,C=matrix(1,3,3)-diag(3),tol=1e-6,scheme="centroid")

plot(testDataEMSB$crit,pch=16,main="RGCCA criterion")
plot(testDataEMSB$stab,pch=16,main="Stability",ylim=c(0,1))
plot(testDataEMSB$obj,pch=16,main="RMSE")

testDataEMSB2=imputeEM(A=A,superblock=TRUE,ncomp=rep(1,3),scale=TRUE,sameBlockWeight=FALSE,tau=rep(1,3),naxis=2,ni=50,C=matrix(1,3,3)-diag(3),tol=1e-6,scheme="centroid")

testDataEMSB3=imputeEM(A=A,superblock=TRUE,ncomp=rep(1,3),scale=TRUE,sameBlockWeight=FALSE,tau=rep(1,3),naxis=3,ni=50,C=matrix(1,3,3)-diag(3),tol=1e-6,scheme="centroid")


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

testDataEM=imputeEM(A=A,ncomp=rep(1,3),superblock=TRUE,scale=TRUE,sameBlockWeight=TRUE,tau=rep(1,3),naxis=1,ni=50,C=matrix(1,3,3)-diag(3),tol=1e-6,scheme="centroid")

plot(testDataEMSB$crit,pch=16)
plot(testDataEMSB$stab,pch=16)

# TO DO test avec deux composantes et un axe
testDataEM=imputeEM(A=A,ncomp=rep(2,3),scale=TRUE,sameBlockWeight=TRUE,tau=rep(1,3),naxis=1,ni=50,C=matrix(1,3,3)-diag(3),tol=1e-8)
testDataEM=imputeEM(A=A,ncomp=rep(2,3),superblock=TRUE,scale=TRUE,sameBlockWeight=TRUE,tau=rep(1,3),naxis=1,ni=50,C=matrix(1,3,3)-diag(3),tol=1e-8)

