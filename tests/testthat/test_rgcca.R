# setwd("/home/caroline.peltier/Bureau/RGCCA")
# lapply(list.files("./R"),function(x){source(paste0("./R/",x))})
data(Russett)
X_agric =as.matrix(Russett[,c("gini","farm","rent")]);
X_ind = as.matrix(Russett[,c("gnpr","labo")]);
X_polit = as.matrix(Russett[ , c("demostab", "dictator")]);
A = list(X_agric);
#C = matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3);

# scaled PCA
resPCA= rgcca (
     blocks=A,
     connection = 1 - diag(length(A)),
     response = NULL,
     superblock = FALSE,
     tau = rep(1, length(A)),
     ncomp = rep(2, length(A)),
     type = "pca",
     verbose = FALSE,
     scheme = "factorial",
     scale = TRUE,
     init = "svd",
     bias = TRUE, 
     tol = 1e-08)

names(resPCA)
(resPCA$astar)
resPCAprcomp=prcomp(A[[1]],scale=TRUE)

varExplPrComp=as.vector((resPCAprcomp$sdev)^2/sum((resPCAprcomp$sdev)^2))[1]
varExplRgcca=resPCA$AVE$AVE_X[[1]][1]
pca_varexpl=round(varExplPrComp-varExplRgcca,digits=4)==0
pca_ind=abs(cor(resPCAprcomp$x[,1],resPCA$Y[[1]][,1]))==1
pca_ind2=abs(cor(resPCAprcomp$x[,2],resPCA$Y[[1]][,2]))==1
pca_var=abs(cor(resPCAprcomp$rotation[,1],resPCA$astar[[1]][,1]))==1
pca_var2=abs(cor(resPCAprcomp$rotation[,2],resPCA$astar[[1]][,2]))==1
test_that("pca_varexpl",{expect_true(pca_varexpl)})
pca_eig=abs(varExplPrComp-varExplRgcca)<1e-8
test_that("pca_ind",{expect_true(pca_ind)})
test_that("pca_var",{expect_true(pca_var)})
test_that("pca_eig",{expect_true(pca_eig)})
test_that("pca_ind2",{expect_true(pca_ind2)})
test_that("pca_var2",{expect_true(pca_var2)})



# unscaled PCA
unscaledPCA= rgcca (
    blocks=A,
    connection = 1 - diag(length(A)),
    response = NULL,
    superblock = FALSE,
    tau = rep(1, length(A)),
    ncomp = rep(2, length(A)),
    type = "pca",
    verbose = FALSE,
    scheme = "factorial",
    scale_block = FALSE,
    scale = FALSE,
    init = "svd",
    bias = TRUE, 
    tol = 1e-08)
unscaledPCAprcomp=prcomp(A[[1]],scale=FALSE)
unscaledvarExplPrComp=as.vector((unscaledPCAprcomp$sdev)^2/sum((unscaledPCAprcomp$sdev)^2))[1]
unscaledvarExplRgcca=unscaledPCA$AVE$AVE_X[[1]][1]
upca_varexpl=round(unscaledvarExplPrComp-unscaledvarExplRgcca,digits=4)==0
upca_ind=abs(cor(unscaledPCAprcomp$x[,1],unscaledPCA$Y[[1]][,1]))==1
upca_var=abs(cor(unscaledPCAprcomp$rotation[,1],unscaledPCA$astar[[1]][,1]))==1
upca_ind2=abs(cor(unscaledPCAprcomp$x[,2],unscaledPCA$Y[[1]][,2]))==1
upca_var2=abs(cor(unscaledPCAprcomp$rotation[,2],unscaledPCA$astar[[1]][,2]))==1

test_that("upca_ind",{expect_true(upca_ind)})
test_that("upca_var",{expect_true(upca_var)})
test_that("upca_ind2",{expect_true(upca_ind)})
test_that("upca_var2",{expect_true(upca_var)})
#test_that("upca_varexpl",{expect_true(upca_varexpl)})

#testthat("upca_eig",{expect_true(abs(unscaledvarExplPrComp-unscaledvarExplRgcca<1e-8))}) #TODO

# With superblock
data(Russett)
X_agric =as.matrix(Russett[,c("gini","farm","rent")]);
X_ind = as.matrix(Russett[,c("gnpr","labo")]);
X_polit = as.matrix(Russett[ , c("demostab", "dictator")]);
A = list(X_agric[,1],X_agric[,2],X_agric[,3]);
scaledPCASB= rgcca (
    blocks=A,
    connection = 1 - diag(length(A)),
    response = NULL,
    superblock = TRUE,
    tau = rep(1, length(A)),
    ncomp = rep(1, length(A)),
    type = "rgcca",
    verbose = TRUE,
    scheme = "factorial",
    scale = TRUE,
    init = "svd",
    bias = TRUE, 
    tol = 1e-08)
pcaSB=prcomp(cbind(X_agric,X_ind,X_polit),scale=TRUE)
pcasb_ind=abs(cor(pcaSB$x[,1],scaledPCASB$Y[[1]][,1]))==1
#pcasb_var=abs(cor(pcaSB$rotation[,1],scaledPCASB$astar[[1]][,1]))==1

#------------PLS 
#  res_pls = plsr(X_polit ~ X_agric, ncomp = 1, method = "simpls")
#  A = list(X_agric,X_polit);
#  pls_with_rgcca= rgcca (
#      blocks=A,
#      connection=matrix(c(0,1,1,0),2,2),
#      tau=rep(1,2),
#      ncomp = rep(1, length(A)),
#     scale_block=FALSE)
#  
# # 
#  cor_X = abs(cor(res_pls$fitted.values[,,1][,1], pls_with_rgcca$Y[[1]]))
#  cor_Y = abs(drop(cor(res_pls$Yloadings, pls_with_rgcca$a[[2]])))
#  cor_X
#  cor_Y
# 
# cor_X
# cor_Y
# library(spls)
# #Loadings are indeed the same (just difference of normalization)
# plot(cor_X*pls_with_rgcca$a[[1]], col = "red", pch = 16, main = "X_loadings", ylab = "loadings")
# points(res_pls$projection/norm2(res_pls$projection))
# 
# plot(cor_X*res_rgcca$a[[2]], col = "red", pch = 16, main = "Y_loadings", ylab = "loadings")
# points(res_pls$Yloadings/norm2(res_pls$Yloadings))

# rgcca (
#     blocks,
#     connection = 1 - diag(length(blocks)),
#     response = NULL,
#     superblock = TRUE,
#     tau = rep(1, length(blocks)),
#     ncomp = rep(2, length(blocks)),
#     type = "rgcca",
#     verbose = TRUE,
#     scheme = "factorial",
#     scale = TRUE,
#     init = "svd",
#     bias = TRUE, 
#     tol = 1e-08)
# 

# Test with block with 1 variable only
 data(Russett)
 X_agric =as.matrix(Russett[,c("gini","farm","rent")]);
 X_ind = as.matrix(Russett[,c("gnpr","labo")]);
 X_polit = as.matrix(Russett[ , c("demostab")]);
 A = list(X_agric,X_ind,X_polit);
 resPCA= rgcca (
     blocks=A,
     connection = 1 - diag(length(A)),
     response = NULL,
     superblock = FALSE,
     tau = rep(1, length(A)),
     ncomp = c(2,2,1),
     type = "rgcca",
     verbose = FALSE,
     scheme = "factorial",
     scale = TRUE,
     init = "svd",
     bias = TRUE, 
     tol = 1e-08)
 
 
 data(Russett)
 X_agric =as.matrix(Russett[,c("gini","farm","rent")]);
 X_ind = as.matrix(Russett[,c("gnpr","labo")]);
 X_polit = as.matrix(Russett[ , c("demostab")]);
 A = list(X_agric,X_ind,X_agric);
 resRGCCA= rgcca (
     blocks=A,
     connection = 1 - diag(length(A)),
     response = NULL,
     superblock = FALSE,
     tau = rep("optimal", length(A)),
     ncomp = rep(2, length(A)),
     type = "rgcca",
     verbose = FALSE,
     scheme = "factorial",
     scale = TRUE,
     init = "svd",
     bias = TRUE, 
     tol = 1e-08)
 
 # Testing quiet=TRUE/quiet=FALSE
 
 data(Russett)
 X_agric =as.matrix(Russett[,c("gini","farm","rent")]);
 X_ind = as.matrix(Russett[,c("gnpr","labo")]);
 X_polit = as.matrix(Russett[ , c("demostab")]);
 A = list(X_agric,X_ind,X_agric);
 names(A)=c("Agri","Ind","Polit")
 resRGCCA= rgcca (
     blocks=A,
     connection = 1 - diag(length(A)),
     response = NULL,
     superblock = FALSE,
     tau = rep("optimal", length(A)),
     ncomp = rep(2, length(A)),
     type = "rgcca",
     verbose = FALSE,
     scheme = "factorial",
     scale = TRUE,
     init = "svd",
     bias = TRUE, 
     tol = 1e-08,quiet=FALSE)
 
 resRGCCA
 
 data(Russett)
 X_agric =as.matrix(Russett[,c("gini","farm","rent")]);
 X_ind = as.matrix(Russett[,c("gnpr","labo")]);
 X_polit = as.matrix(Russett[ , 6:11]);
 A = list(X_agric,X_ind,X_polit);
 names(A)=c("Agri","Ind","Polit")
 
 C0=matrix(0,3,3);C0[2:3,1]=1;C0[1,2:3]=1
 C1=matrix(0,3,3);C1[1:2,3]=1;C1[3,1:2]=1
 A1=list(A[[2]],A[[3]],A[[1]])
 resRgccaNipals3=rgcca(blocks=A1,connection=C1,type="rgcca",method="nipals",ncomp=2)
 resRgccaNipals=rgcca(blocks=A,connection=C0,type="rgcca",method="nipals",ncomp=2)
 head(resRgccaNipals3$Y[[3]])
 head(resRgccaNipals$Y[[1]])
 
 
 mat=matrix(rnorm(500),nrow=10,ncol=50);rownames(mat)=paste0("S",1:10);colnames(mat)=paste0("R",1:50)
 A=list(mat)
 resPCA= rgcca (
     blocks=A,
     connection = 1 - diag(length(A)),
     response = NULL,
     superblock = FALSE,
     tau = rep(1, length(A)),
     ncomp = rep(2, length(A)),
     type = "pca",
     verbose = FALSE,
     scheme = "factorial",
     scale = TRUE,
     init = "svd",
     bias = TRUE, 
     tol = 1e-08)
 
 names(resPCA)
 (resPCA$astar)
 resPCAprcomp=prcomp(A[[1]],scale=TRUE)
 pca_ind=abs(cor(resPCAprcomp$x[,1],resPCA$Y[[1]][,1]))==1
 
 