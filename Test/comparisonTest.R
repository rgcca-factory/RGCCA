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
A = list(agri=X_agric, ind=X_ind, polit=X_polit)

refRgcca=rgcca(Aref, C = 1 - diag(length(A)), tau = rep(1, length(A)), ncomp = rep(2, length(A)), scheme = "centroid", scale = TRUE,   init = "svd", bias = TRUE, tol = 1e-08, verbose = TRUE,sameBlockWeight=TRUE,returnA=TRUE,na.rm=FALSE) 
listNAdataset=  createNA(A=Aref,option="block",pNA=rep(0.2,3),nAllRespondants=10,output="list")
selectCompletePatient=listNAdataset$subjectKept
methodRgcca=rgccaNa(A=listNAdataset$dat,method="em",verbose=FALSE,ncomp=rep(2,length(A)),returnA=TRUE)
comparison(rgcca1=refRgcca,rgcca2=methodRgcca$rgcca,selectPatient=selectCompletePatient,indNA=methodRgcca$indNA)
