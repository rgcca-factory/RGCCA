# '# ImputeSB test
# 
# '''
#setwd("/home/caroline.peltier/Bureau/RGCCA")
rm(list=ls())
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
testDataSB=imputeSB(A=A,ncomp=rep(1,3),scale=TRUE,sameBlockWeight=TRUE,tau=rep(1,3),naxis=1,ni=4)
plot(testDataSB$crit,pch=16)
plot(testDataSB$stab,pch=16)


source("/home/caroline.peltier/Bureau/EtudeNA/readDataset.r")
repertoire=paste("/home/caroline.peltier/Bureau/EtudeNA/Datasets/Biosca/bloc5/",sep="")
setwd("/home/caroline.peltier/Bureau/RGCCA")
namesFiles=dir("./R")
# loading functions in R directory
sapply(namesFiles,function(x){source(paste0("./R/",x))})

j=1
nAxe=2;nBlock=3;scale=TRUE;sameBlockWeight=TRUE,tau=rep(1,3)
 setwd(paste(repertoire,j,sep=""))
  testData=readDataset(c("CLI","MRS","VOL"))
 resimpute= imputeSB(testData,ncomp=rep(nAxe,nBlock),scale=scale,sameBlockWeight=sameBlockWeight,tau=rep(1,3),tol=1e-8,ni=50,naxis=1)
  plot(resimpute$crit)
  plot(resimpute$obj)