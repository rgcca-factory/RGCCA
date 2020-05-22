# '# Test rgcca
# 
# '''
#setwd("~/Bureau/RGCCA/tests/testthat")
data(Russett)
X_agric =as.matrix(Russett[,c("gini","farm","rent")]);
X_ind = as.matrix(Russett[,c("gnpr","labo")]);
X_polit = as.matrix(Russett[ , c("demostab", "dictator")]);
A = list(X_agric, X_ind, X_polit);
C = matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3);

resSgcca = sgcca(A, C, ncomp=rep(2,3),sparsity= c(1, 1, 1), scheme = "factorial", scale = TRUE,verbose=FALSE)

resSgcca=sgcca(A,C,sparsity=rep(0.8,3))
head(resSgcca$Y[[3]])
A1=list(pol=A[[3]],agr=A[[1]],ind=A[[2]])
C1=matrix(c(0,1,1,1,0,0,1,0,0),3,3)
resSgcca1=sgcca(A1,C1,sparsity=rep(0.8,3))
head(resSgcca1$Y[[1]])


resRgcca=rgcca(blocks=A,connection=C)
head(resRgcca$Y[[3]])
resRgcca1=rgcca(blocks=A1,connection=C1)
head(resRgcca1$Y[[1]])
