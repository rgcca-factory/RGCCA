# '# writeList test
# 
# '''
#setwd("/home/caroline.peltier/Bureau/RGCCA")
data(Russett)
library(missMDA)
# Russetts
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
res1=rgccaNa (A,method="sem", C = 1 - diag(length(A)), tau = rep(1, length(A)), scale=TRUE,  ncomp = rep(1, length(A)), scheme = "centroid",   init = "svd", bias = TRUE, tol = 1e-08, verbose = FALSE,sameBlockWeight=TRUE,returnA=FALSE,knn.k="all",knn.output="weightedMean",knn.klim=NULL,knn.sameBlockWeight=TRUE) 
res2=rgccaNa (A,method="mean", C = 1 - diag(length(A)), tau = rep(1, length(A)),  ncomp = rep(1, length(A)), scheme = "centroid", scale = TRUE,   init = "svd", bias = TRUE, tol = 1e-08, verbose = FALSE,sameBlockWeight=TRUE,returnA=FALSE,knn.k="all",knn.output="weightedMean",knn.klim=NULL,knn.sameBlockWeight=TRUE) 
res3=rgccaNa (A,method="em", C = 1 - diag(length(A)), tau = rep(1, length(A)),   ncomp = rep(1, length(A)), scheme = "centroid", scale = TRUE,   init = "svd", bias = TRUE, tol = 1e-05, verbose = FALSE,sameBlockWeight=TRUE,returnA=FALSE,knn.k="all",knn.output="weightedMean",knn.klim=NULL,knn.sameBlockWeight=TRUE) 
res4=rgccaNa (A,method="sem3", C = 1 - diag(length(A)), tau = rep(1, length(A)),    ncomp = rep(1, length(A)), scheme = "centroid", scale = TRUE,   init = "svd", bias = TRUE, tol = 1e-08, verbose = FALSE,sameBlockWeight=TRUE,returnA=FALSE,knn.k="all",knn.output="weightedMean",knn.klim=NULL,knn.sameBlockWeight=TRUE) 
res5=rgccaNa (A,method="pca", C = 1 - diag(length(A)), tau = rep(1, length(A)),     ncomp = rep(1, length(A)), scheme = "centroid", scale = TRUE,   init = "svd", bias = TRUE, tol = 1e-08, verbose = FALSE,sameBlockWeight=TRUE,returnA=FALSE,knn.k="all",knn.output="weightedMean",knn.klim=NULL,knn.sameBlockWeight=TRUE) 
res6=rgccaNa (A,method="rpca", C = 1 - diag(length(A)), tau = rep(1, length(A)),    ncomp = rep(1, length(A)), scheme = "centroid", scale = TRUE,   init = "svd", bias = TRUE, tol = 1e-08, verbose = FALSE,sameBlockWeight=TRUE,returnA=FALSE,knn.k="all",knn.output="weightedMean",knn.klim=NULL,knn.sameBlockWeight=TRUE) 
res7=rgccaNa (A,method="mfa", C = 1 - diag(length(A)), tau = rep(1, length(A)),  ncomp = rep(1, length(A)), scheme = "centroid", scale = TRUE,   init = "svd", bias = TRUE, tol = 1e-08, verbose = FALSE,sameBlockWeight=TRUE,returnA=FALSE,knn.k="all",knn.output="weightedMean",knn.klim=NULL,knn.sameBlockWeight=TRUE) 
#res8=rgccaNa (A,method="imputeInRgcca", C = 1 - diag(length(A)), tau = rep(1, length(A)),    ncomp = rep(1, length(A)), scheme = "centroid", scale = TRUE,   init = "svd", bias = TRUE, tol = 1e-08, verbose = FALSE,sameBlockWeight=TRUE,returnA=FALSE,knn.k="all",knn.output="weightedMean",knn.klim=NULL,knn.sameBlockWeight=TRUE) 
res9=rgccaNa (A,method="imputeInRgccaLL", C = 1 - diag(length(A)), tau = rep(1, length(A)),    ncomp = rep(1, length(A)), scheme = "centroid", scale = TRUE,   init = "svd", bias = TRUE, tol = 1e-08, verbose = FALSE,sameBlockWeight=TRUE,returnA=FALSE,knn.k="all",knn.output="weightedMean",knn.klim=NULL,knn.sameBlockWeight=TRUE) 
 test_that("rgccaNa_1",{expect_true(is.list(res1))})
 test_that("rgccaNa_2",{expect_true(is.list(res2))})
 test_that("rgccaNa_3",{expect_true(is.list(res3))})
 test_that("rgccaNa_4",{expect_true(is.list(res4))})
 test_that("rgccaNa_5",{expect_true(is.list(res5))})
 test_that("rgccaNa_6",{expect_true(is.list(res6))})
 test_that("rgccaNa_7",{expect_true(is.list(res7))})
 test_that("rgccaNa_9",{expect_true(is.list(res9))})
 
