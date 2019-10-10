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
testDataSB=imputeSB(A=A,ncomp=rep(1,3),scale=TRUE,sameBlockWeight=TRUE,tau=rep(1,3),naxis=1,ni=4)
plot(testDataSB$crit,pch=16)
plot(testDataSB$stab,pch=16)
