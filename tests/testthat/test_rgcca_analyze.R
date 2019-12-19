# setwd("/home/caroline.peltier/Bureau/RGCCA")
# lapply(list.files("./R"),function(x){source(paste0("./R/",x))})
data(Russett)
X_agric =as.matrix(Russett[,c("gini","farm","rent")]);
X_ind = as.matrix(Russett[,c("gnpr","labo")]);
X_polit = as.matrix(Russett[ , c("demostab", "dictator")]);
A = list(X_agric);
#C = matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3);

# scaled PCA
resPCA= rgcca.analyze (
     blocks=A,
     connection = 1 - diag(length(A)),
     response = NULL,
     superblock = TRUE,
     tau = rep(1, length(A)),
     ncomp = rep(2, length(A)),
     type = "pca",
     verbose = TRUE,
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

test_that("pca_ind",{expect_true(cor(resPCAprcomp$x[,1],resPCA$Y[[1]][,1])==1)})
test_that("pca_var",{expect_true(cor(resPCAprcomp$rotation[,1],resPCA$astar[[1]][,1])==1)})
#test_that("pca_eig",{expect_true(abs(varExplPrComp-varExplRgcca<1e-8))})

# unscaled PCA
unscaledPCA= rgcca.analyze (
    blocks=A,
    connection = 1 - diag(length(A)),
    response = NULL,
    superblock = TRUE,
    tau = rep(1, length(A)),
    ncomp = rep(2, length(A)),
    type = "pca",
    verbose = TRUE,
    scheme = "factorial",
    scale = FALSE,
    init = "svd",
    bias = TRUE, 
    tol = 1e-08)
unscaledPCAprcomp=prcomp(A[[1]],scale=FALSE)
unscaledvarExplPrComp=as.vector((unscaledPCAprcomp$sdev)^2/sum((unscaledPCAprcomp$sdev)^2))[1]
unscaledvarExplRgcca=unscaledPCA$AVE$AVE_X[[1]][1]

test_that("upca_ind",{expect_true(cor(unscaledPCAprcomp$x[,1],unscaledPCA$Y[[1]][,1])==1)})
test_that("upca_var",{expect_true(cor(unscaledPCAprcomp$rotation[,1],unscaledPCA$astar[[1]][,1])==1)})
#testthat("upca_eig",{expect_true(abs(unscaledvarExplPrComp-unscaledvarExplRgcca<1e-8))})
#------------PLS 
# library(pls)
# res_pls = plsr(X_polit ~ X_agric, ncomp = 1, method = "simpls")
# A = list(X_agric,X_polit);
# pls_with_rgcca= rgcca.analyze (
#     blocks=A,
#     ncomp = rep(1, length(A)),
#     type = "pls",sameBlockWeight=TRUE)
# 
# 
# cor_X = drop(cor(res_pls$projection, pls_with_rgcca$a[[1]]))
# cor_Y = drop(cor(res_pls$Yloadings, pls_with_rgcca$a[[2]]))
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

# rgcca.analyze (
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

# 
