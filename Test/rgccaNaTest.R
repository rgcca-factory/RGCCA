# '# writeList test
# 
# '''
#setwd("/home/caroline.peltier/Bureau/RGCCA")
namesFiles=dir("./R")
sapply(namesFiles,function(x){source(paste0("./R/",x))})

# Russetts
library(missMDA)
X_agric =as.matrix(Russett[,c("gini","farm","rent")])
X_ind = as.matrix(Russett[,c("gnpr","labo")])
X_polit = as.matrix(Russett[ , c("demostab", "dictator")])
X_agric[c(2,4),]=NA
X_ind[1,]=NA
X_polit[5,1]=NA
A = list(agri=X_agric, ind=X_ind, polit=X_polit)
A_ref=list(agri=as.matrix(Russett[,c("gini","farm","rent")]),ind=as.matrix(Russett[,c("gnpr","labo")]),polit=as.matrix(Russett[ , c("demostab", "dictator")]))
A_ref2=lapply(A_ref,scale)
A2=lapply(A,scale)

rgccaNa (A,method="mean", C = 1 - diag(length(A)), tau = rep(1, length(A)), refData=NULL,    ncomp = rep(1, length(A)), scheme = "centroid", scale = TRUE,   init = "svd", bias = TRUE, tol = 1e-08, verbose = TRUE,na.impute="none",na.niter=10,na.keep=NULL,nboot=10,sameBlockWeight=TRUE,returnA=FALSE,knn.k="all",knn.output="weightedMean",knn.klim=NULL,knn.sameBlockWeight=TRUE) 
rgccaNa (A,method="em", C = 1 - diag(length(A)), tau = rep(1, length(A)), refData=NULL,    ncomp = rep(1, length(A)), scheme = "centroid", scale = TRUE,   init = "svd", bias = TRUE, tol = 1e-08, verbose = TRUE,na.impute="none",na.niter=10,na.keep=NULL,nboot=10,sameBlockWeight=TRUE,returnA=FALSE,knn.k="all",knn.output="weightedMean",knn.klim=NULL,knn.sameBlockWeight=TRUE) 
rgccaNa (A,method="sem3", C = 1 - diag(length(A)), tau = rep(1, length(A)), refData=NULL,    ncomp = rep(1, length(A)), scheme = "centroid", scale = TRUE,   init = "svd", bias = TRUE, tol = 1e-08, verbose = TRUE,na.impute="none",na.niter=10,na.keep=NULL,nboot=10,sameBlockWeight=TRUE,returnA=FALSE,knn.k="all",knn.output="weightedMean",knn.klim=NULL,knn.sameBlockWeight=TRUE) 
rgccaNa (A,method="pca", C = 1 - diag(length(A)), tau = rep(1, length(A)), refData=NULL,    ncomp = rep(1, length(A)), scheme = "centroid", scale = TRUE,   init = "svd", bias = TRUE, tol = 1e-08, verbose = TRUE,na.impute="none",na.niter=10,na.keep=NULL,nboot=10,sameBlockWeight=TRUE,returnA=FALSE,knn.k="all",knn.output="weightedMean",knn.klim=NULL,knn.sameBlockWeight=TRUE) 
rgccaNa (A,method="rpca", C = 1 - diag(length(A)), tau = rep(1, length(A)), refData=NULL,    ncomp = rep(1, length(A)), scheme = "centroid", scale = TRUE,   init = "svd", bias = TRUE, tol = 1e-08, verbose = TRUE,na.impute="none",na.niter=10,na.keep=NULL,nboot=10,sameBlockWeight=TRUE,returnA=FALSE,knn.k="all",knn.output="weightedMean",knn.klim=NULL,knn.sameBlockWeight=TRUE) 
rgccaNa (A,method="mfa", C = 1 - diag(length(A)), tau = rep(1, length(A)), refData=NULL,    ncomp = rep(1, length(A)), scheme = "centroid", scale = TRUE,   init = "svd", bias = TRUE, tol = 1e-08, verbose = TRUE,na.impute="none",na.niter=10,na.keep=NULL,nboot=10,sameBlockWeight=TRUE,returnA=FALSE,knn.k="all",knn.output="weightedMean",knn.klim=NULL,knn.sameBlockWeight=TRUE) 
rgccaNa (A,method="imputeInRgcca", C = 1 - diag(length(A)), tau = rep(1, length(A)), refData=NULL,    ncomp = rep(1, length(A)), scheme = "centroid", scale = TRUE,   init = "svd", bias = TRUE, tol = 1e-08, verbose = TRUE,na.impute="none",na.niter=10,na.keep=NULL,nboot=10,sameBlockWeight=TRUE,returnA=FALSE,knn.k="all",knn.output="weightedMean",knn.klim=NULL,knn.sameBlockWeight=TRUE) 
