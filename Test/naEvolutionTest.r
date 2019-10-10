setwd("/home/caroline.peltier/Bureau/RGCCA")
rm(list=ls())
library(RGCCA)
library(MASS)
library(nipals)
namesFiles=dir("./R")
# loading functions in R directory
sapply(namesFiles,function(x){source(paste0("./R/",x))})
library(parallel)
set.seed(42);X1=matrix(rnorm(350),70,5);X2=matrix(rnorm(280),70,4)
A=list(X1,X2)
listResults=naEvolution(A=A,prctNA=c(0.1,0.2,0.3,0.4),listMethods=c("mean","complete"))
plotEvol(listResults,ylim=NULL)

#--------------------
# test on russets
#---------------------
X_agric =as.matrix(Russett[,c("gini","farm","rent")])
X_ind = as.matrix(Russett[,c("gnpr","labo")])
X_polit = as.matrix(Russett[ , c("demostab", "dictator")])
X_agric[c(2,4),]=NA
X_ind[1,]=NA
X_polit[5,1]=NA
A = list(agri=X_agric, ind=X_ind, polit=X_polit)
listResults=naEvolution(A=A,listMethods=c("complete","mean","em","superblockEM","knn4","nipals"),prctNA=c(0.1,0.2,0.3))
# bug pour plus grand que 04
plotEvol(listResults,ylim=c(0.9,1),output="a")
plotEvol(listResults,ylim=NULL,output="rv")
plotEvol(listResults,ylim=NULL,output="rvComplete")
plotEvol(listResults,ylim=NULL,output="bm")

# Tests for Congress PLS
data(Russett)
library(missMDA)
library(FactoMineR)
library(parallel)
X_agric =as.matrix(Russett[,c("gini","farm","rent")])
X_ind = as.matrix(Russett[,c("gnpr","labo")])
X_polit = as.matrix(Russett[ , c("demostab", "dictator")])
A = list(agri=X_agric, ind=X_ind, polit=X_polit)
#ponctual
listResults=naEvolution(A=A,listMethods=c("complete","nipals","imputeInRgcca"),prctNA=c(0.05,0.1,0.15,0.2,0.25,0.3),typeNA="ponc",ncomp=rep(1,3))
plotEvol(listResults,output="a",barType = "stderr",ylim=c(0,0.3))
plotEvol(listResults,ylim=c(0.6,1),output="rv")
plotEvol(listResults,ylim=c(0.9,1),output="rvComplete",barType = "stderr")
# block
listResults=naEvolution(A=A,listMethods=c("complete","nipals","imputeInRgcca"),prctNA=c(0.05,0.1,0.15,0.2,0.25,0.3),typeNA="block",ncomp=rep(1,3))
plotEvol(listResults,ylim=c(0.9,1),output="a",block=1,barType = "stderr")
plotEvol(listResults,ylim=c(0.6,1),output="rv")
plotEvol(listResults,ylim=c(0.9,1),output="rvComplete",barType = "stderr")

# differences ponctuelles
listResults=naEvolution(A=A,listMethods=c("pca","complete","mean","nipals","iterativeSB","em","sem","imputeInRgcca"),prctNA=c(0.1,0.2,0.3),typeNA="ponc")
plotEvol(listResults,ylim=c(0.0,1),output="a",block=1,barType = "stderr")
plotEvol(listResults,ylim=c(0.6,1),output="rv")
plotEvol(listResults,ylim=c(0.95,1),output="rvComplete")
plotEvol(listResults,ylim=c(0,1),output="rmse")

#------------------
# test on biosca
#------------------
setwd("/home/caroline.peltier/Bureau/EtudeNA/Datasets/Biosca/Reference")
refData=readDataset(c("CLI","MRS","VOL"))
listResults=naEvolution(A=refData,listMethods=c("complete","mean","em","sem1","knn4","nipals"),prctNA=c(0.1,0.2,0.3))
A=refData

A2=imputeEM(A=A,superblock=TRUE,ncomp=ncomp,scale=scale,sameBlockWeight=sameBlockWeight,tau=tau,naxis=as.numeric(substr(method,4,4)),ni=50,C=C,tol=1e-6)$A
plotEvol(listResults,ylim=NULL,output="a")
plotEvol(listResults,ylim=NULL,output="rv")
plotEvol(listResults,ylim=NULL,output="rvComplete")
plotEvol(listResults,ylim=NULL,output="bm")
plotEvol(listResults,ylim=c(0,1),output="rmse")
