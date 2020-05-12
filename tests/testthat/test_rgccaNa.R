# # '# writeList test
# # 
# # '''
# #setwd("/home/caroline.peltier/Bureau/RGCCA")
 data(Russett)
# #library(missMDA)
# # Russetts
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
 res1=rgccaNa (A,method="sem", connection = 1 - diag(length(A)), tau = rep(1, length(A)), scale=TRUE,  ncomp = rep(1, length(A)), scheme = "centroid",   init = "svd", bias = TRUE, tol = 1e-08, verbose = FALSE,scale_block=TRUE,knn.k="all",knn.output="weightedMean",knn.klim=NULL,knn.scale_block=TRUE)
 res2=rgccaNa (A,method="mean", connection = 1 - diag(length(A)), tau = rep(1, length(A)),  ncomp = rep(1, length(A)), scheme = "centroid", scale = TRUE,   init = "svd", bias = TRUE, tol = 1e-08, verbose = FALSE,scale_block=TRUE,knn.k="all",knn.output="weightedMean",knn.klim=NULL,knn.scale_block=TRUE)
 res3=rgccaNa (A,method="em", connection = 1 - diag(length(A)), tau = rep(1, length(A)),   ncomp = rep(1, length(A)), scheme = "centroid", scale = TRUE,   init = "svd", bias = TRUE, tol = 1e-05, verbose = FALSE,scale_block=TRUE,knn.k="all",knn.output="weightedMean",knn.klim=NULL,knn.scale_block=TRUE)
 res4=rgccaNa (A,method="sem3", connection = 1 - diag(length(A)), tau = rep(1, length(A)),    ncomp = rep(1, length(A)), scheme = "centroid", scale = TRUE,   init = "svd", bias = TRUE, tol = 1e-08, verbose = FALSE,scale_block=TRUE,knn.k="all",knn.output="weightedMean",knn.klim=NULL,knn.scale_block=TRUE)
 #res5=rgccaNa (A,method="pca", C = 1 - diag(length(A)), tau = rep(1, length(A)),     ncomp = rep(1, length(A)), scheme = "centroid", scale = TRUE,   init = "svd", bias = TRUE, tol = 1e-08, verbose = FALSE,scale_block=TRUE,knn.k="all",knn.output="weightedMean",knn.klim=NULL,knn.scale_block=TRUE)
 #res6=rgccaNa (A,method="rpca", C = 1 - diag(length(A)), tau = rep(1, length(A)),    ncomp = rep(1, length(A)), scheme = "centroid", scale = TRUE,   init = "svd", bias = TRUE, tol = 1e-08, verbose = FALSE,scale_block=TRUE,knn.k="all",knn.output="weightedMean",knn.klim=NULL,knn.scale_block=TRUE)
 #res7=rgccaNa (A,method="mfa", C = 1 - diag(length(A)), tau = rep(1, length(A)),  ncomp = rep(1, length(A)), scheme = "centroid", scale = TRUE,   init = "svd", bias = TRUE, tol = 1e-08, verbose = FALSE,scale_block=TRUE,knn.k="all",knn.output="weightedMean",knn.klim=NULL,knn.scale_block=TRUE)
 #res8=rgccaNa (A,method="imputeInRgcca", C = 1 - diag(length(A)), tau = rep(1, length(A)),    ncomp = rep(1, length(A)), scheme = "centroid", scale = TRUE,   init = "svd", bias = TRUE, tol = 1e-08, verbose = FALSE,scale_block=TRUE,knn.k="all",knn.output="weightedMean",knn.klim=NULL,knn.scale_block=TRUE)
 res9=rgccaNa (A,method="imputeInRgccaLL", connection = 1 - diag(length(A)), tau = rep(1, length(A)),    ncomp = rep(1, length(A)), scheme = "centroid", scale = TRUE,   init = "svd", bias = TRUE, tol = 1e-08, verbose = FALSE,scale_block=TRUE,knn.k="all",knn.output="weightedMean",knn.klim=NULL,knn.scale_block=TRUE)
  test_that("rgccaNa_1",{expect_true(is.list(res1))})
  test_that("rgccaNa_2",{expect_true(is.list(res2))})
  test_that("rgccaNa_3",{expect_true(is.list(res3))})
  test_that("rgccaNa_4",{expect_true(is.list(res4))})
 # test_that("rgccaNa_5",{expect_true(is.list(res5))})
  #test_that("rgccaNa_6",{expect_true(is.list(res6))})
  #test_that("rgccaNa_7",{expect_true(is.list(res7))})
  test_that("rgccaNa_9",{expect_true(is.list(res9))})
 

  X1=matrix(rnorm(150),10,15)
  X2=matrix(rnorm(150),10,15)
  A = list(X1=X1,X2=X2)
  res1=rgccaNa (A,method="sem", connection = 1 - diag(length(A)), tau = rep(1, length(A)), scale=TRUE,  ncomp = rep(1, length(A)), scheme = "centroid",   init = "svd", bias = TRUE, tol = 1e-08, verbose = FALSE,scale_block=TRUE,knn.k="all",knn.output="weightedMean",knn.klim=NULL,knn.scale_block=TRUE)
  res2=rgccaNa (A,method="em", connection = 1 - diag(length(A)), tau = rep(1, length(A)),  ncomp = rep(1, length(A)), scheme = "centroid", scale = TRUE,   init = "svd", bias = TRUE, tol = 1e-08, verbose = FALSE,scale_block=TRUE,knn.k="all",knn.output="weightedMean",knn.klim=NULL,knn.scale_block=TRUE)
  