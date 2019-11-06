# setwd("/home/caroline.peltier/Bureau/RGCCA")
rm(list=ls())
library(RGCCA)
library(MASS)
library(nipals)
library(parallel)
namesFiles=dir("./R")
namesFiles2=namesFiles[!namesFiles%in%c("find.biomarkers.R","optspars_cv.R","plotOptpars3D.R","plotOptspars.R","weight.bootstrap.R")]
# loading functions in R directory
#sapply(namesFiles,function(x){source(paste0("./R/",x))})

for(i in 1:length(namesFiles2))
{
    print(namesFiles2[i])
    source(paste0("./R/",namesFiles2[i]))
}

set.seed(42);X1=matrix(rnorm(350),70,5);X2=matrix(rnorm(280),70,4)
A=list(X1,X2)
listResults=naEvolution(refData=A,prctNA=c(0.1,0.2,0.3,0.4),listMethods=c("mean","complete","nipals","em","sem1","knn4","complete"))
plotEvol(listFinale=listResults,ylim=c(0,1),output="a")

#--------------------
# test on russets
#---------------------
library(MASS)
library(FactoMineR)
library(parallel)
library(missMDA)
data(Russett)
X_agric =as.matrix(Russett[,c("gini","farm","rent")])
X_ind = as.matrix(Russett[,c("gnpr","labo")])
X_polit = as.matrix(Russett[ , colnames(Russett)%in%c("demostab", "dictatur","dictator")])
A = list(agri=X_agric, ind=X_ind, polit=X_polit)

# TODO
# 1- Corriger le superbloc en scale = FALSE et sameBlockWeight=TRUE
# 2- 
# 3- 
# pour le cas le plus simple : tau=1 scale = TRUE, sameBlockWeight=TRUE
listResults=naEvolution(A=A,listMethods=c("complete","nipals","em","pca","sem","imputeInRgccaLL"),prctNA=c(0.1,0.2,0.3,0.4),typeNA = "ponc",nDatasets=10,sameBlockWeight = TRUE,scale=TRUE,tol=1e-6,verbose=TRUE,scheme="horst")

listResults=naEvolution(A=A,listMethods=c("complete","nipals","em","sem","pca"),prctNA=c(0.05,0.1,0.15,0.2,0.25),typeNA = "ponc",nDatasets=20,sameBlockWeight = TRUE,scale=TRUE,tol=1e-6,verbose=TRUE,scheme="horst")

plotEvol(listResults,ylim=c(0,0.3),output="a")

listResults=naEvolution(A=A,listMethods=c("complete","nipals","em","pca","sem"),prctNA=c(0.1,0.2,0.3),typeNA = "ponc",nDatasets=1,sameBlockWeight = TRUE,scale=TRUE,tol=1e-6,verbose=TRUE)


plotEvol(listResults,ylim=c(0,2),output="rmse")
# le cas EM tourne mais et donne des résultats corrects!

#  tau=1 scale = TRUE, sameBlockWeight=FALSE
listResults=naEvolution(A=A,listMethods=c("complete","nipals","em","pca"),prctNA=c(0.1,0.2),typeNA = "ponc",nDatasets=1,sameBlockWeight = FALSE,scale=TRUE,tol=1e-6)
plotEvol(listResults,ylim=c(0,0.2),output="a")
# tout donne de mauvais résultats ??? 

# pour le cas le plus simple : tau=1 scale = FALSE, sameBlockWeight=FALSE
listResults=naEvolution(A=A,listMethods=c("complete","nipals","em","pca"),prctNA=c(0.1,0.2),typeNA = "ponc",nDatasets=1,sameBlockWeight = FALSE,scale=FALSE,tol=1e-6)
plotEvol(listResults,ylim=c(0,0.5),output="a")
# donne de mauvais résultats ! ! !

# pour le cas le plus simple : tau=1 scale = FALSE, sameBlockWeight=TRUE
listResults=naEvolution(A=A,listMethods=c("complete","nipals","em","pca"),prctNA=c(0.1,0.2,0.3),typeNA = "ponc",nDatasets=1,sameBlockWeight = TRUE,scale=FALSE,tol=1e-6)
plotEvol(listResults,ylim=c(0,0.2),output="a")
# le em  fonctionne bien !




listResults=naEvolution(A=A,listMethods=c("complete","em","sem","pca"),prctNA=c(0.1,0.2,0.3),typeNA = "ponc",nDatasets=1,sameBlockWeight = TRUE,scale=FALSE,tol=1e-6)

# pour le scale = FALSE, sameBlockWeight=TRUE : NON PROGRAMME POUR EM
listResults=naEvolution(A=A,listMethods=c("complete","nipals","pca"),prctNA=c(0.1,0.2,0.3),typeNA = "ponc",nDatasets=1,sameBlockWeight = TRUE,scale=FALSE,tol=1e-6)

listResults=naEvolution(A=A,listMethods=c("complete","em","sem","pca"),prctNA=c(0.1,0.2,0.3),typeNA = "ponc",nDatasets=1,sameBlockWeight = FALSE,scale=TRUE,tol=1e-6,verbose=TRUE)
# bug pour plus grand que 04
plotEvol(listResults,ylim=c(0,0.12),output="a")

plotEvol(listResults,ylim=NULL,output="rv")
plotEvol(listResults,ylim=NULL,output="rvComplete")
plotEvol(listResults,ylim=NULL,output="bm")
#------------------
# test on biosca
#------------------
setwd("/home/caroline.peltier/Bureau/EtudeNA/Datasets/Biosca/Reference")
refData=readDataset(c("CLI","MRS","VOL"))
listResults=naEvolution(A=refData,listMethods=c("complete","mean","em","sem1","knn4","nipals"),prctNA=c(0.1,0.2,0.3),typeNA="ponc")
A=refData
plotEvol(listResults,ylim=NULL,output="a")
plotEvol(listResults,ylim=NULL,output="rv")
plotEvol(listResults,ylim=NULL,output="rvComplete")
plotEvol(listResults,ylim=NULL,output="bm")
plotEvol(listResults,ylim=c(0.4,1.5),output="rmse")
