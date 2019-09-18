# '# writeList test
# 
# '''
set.seed(42);X1=matrix(rnorm(350),70,5);X2=matrix(rnorm(280),70,4)
A=list(X1,X2)

rgccaNa (A,method="mean", C = 1 - diag(length(A)), tau = rep(1, length(A)), refData=NULL,    ncomp = rep(1, length(A)), scheme = "centroid", scale = TRUE,   init = "svd", bias = TRUE, tol = 1e-08, verbose = TRUE,na.impute="none",na.niter=10,na.keep=NULL,nboot=10,sameBlockWeight=TRUE,returnA=FALSE,knn.k="all",knn.output="weightedMean",knn.klim=NULL,knn.sameBlockWeight=TRUE) 

rgccaNa (A,method="em", C = 1 - diag(length(A)), tau = rep(1, length(A)), refData=NULL,    ncomp = rep(1, length(A)), scheme = "centroid", scale = TRUE,   init = "svd", bias = TRUE, tol = 1e-08, verbose = TRUE,na.impute="none",na.niter=10,na.keep=NULL,nboot=10,sameBlockWeight=TRUE,returnA=FALSE,knn.k="all",knn.output="weightedMean",knn.klim=NULL,knn.sameBlockWeight=TRUE) 
rgccaNa (A,method="superblockEM", C = 1 - diag(length(A)), tau = rep(1, length(A)), refData=NULL,    ncomp = rep(1, length(A)), scheme = "centroid", scale = TRUE,   init = "svd", bias = TRUE, tol = 1e-08, verbose = TRUE,na.impute="none",na.niter=10,na.keep=NULL,nboot=10,sameBlockWeight=TRUE,returnA=FALSE,knn.k="all",knn.output="weightedMean",knn.klim=NULL,knn.sameBlockWeight=TRUE) 
