# '# Test rgccak
# 
# '''
# setwd("./tests/testthat")
library(MASS)
data(Russett)
X_agric =as.matrix(Russett[,c("gini","farm","rent")]);
X_ind = as.matrix(Russett[,c("gnpr","labo")]);
X_polit = as.matrix(Russett[ , c("demostab", "dictator")]);
A = list(X_agric, X_ind, X_polit);
C = matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3);

resultRgccak_Tau1_test = rgccak(A, C, tau = c(1, 1, 1), scheme = "factorial")
prodYA1_test=resultRgccak_Tau1_test$Y[,1]%*%t(resultRgccak_Tau1_test$a[[1]])
prodYA2_test=resultRgccak_Tau1_test$Y[,2]%*%t(resultRgccak_Tau1_test$a[[2]])
prodYA3_test=resultRgccak_Tau1_test$Y[,3]%*%t(resultRgccak_Tau1_test$a[[3]])

# test_that("test_rgcca_tau1",{
#     load(file="../results/resultRgccak_Tau1")
#     expect_true(all.equal(resultRgccak_Tau1_test[3:7],resultRgccak_Tau1[3:7]))
# })

test_that("test_rgcca_tau1ya1",{
    load(file="../results/resultRgccak_Tau1")
    prodYA1=resultRgccak_Tau1$Y[,1]%*%t(resultRgccak_Tau1$a[[1]])
     expect_true(all.equal(prodYA1,prodYA1_test))
  })

test_that("test_rgcca_tau1ya2",{
    load(file="../results/resultRgccak_Tau1")
    prodYA2=resultRgccak_Tau1$Y[,2]%*%t(resultRgccak_Tau1$a[[2]])
    expect_true(all.equal(prodYA2,prodYA2_test))
})
test_that("test_rgcca_tau1ya3",{
    load(file="../results/resultRgccak_Tau1")
    prodYA3=resultRgccak_Tau1$Y[,3]%*%t(resultRgccak_Tau1$a[[3]])
    expect_true(all.equal(prodYA3,prodYA3_test))
})
# Tau 0
#------------------
load(file="../results/resultRgccak_Tau0")
resultRgccak_Tau0_test = rgccak(A, C, tau = c(0, 0, 0), scheme = "factorial")

# test_that("test_rgccak_tau0",{
#     load(file="../results/resultRgccak_Tau0")
#    # testing the same results than previously
#     expect_true(all.equal(resultRgccak_Tau0_test[3:7],resultRgccak_Tau0[3:7]))
# })
prodYA1_test=resultRgccak_Tau0_test$Y[,1]%*%t(resultRgccak_Tau0_test$a[[1]])
prodYA2_test=resultRgccak_Tau0_test$Y[,2]%*%t(resultRgccak_Tau0_test$a[[2]])
prodYA3_test=resultRgccak_Tau0_test$Y[,3]%*%t(resultRgccak_Tau0_test$a[[3]])

test_that("test_rgcca_tau0ya1",{
    load(file="../results/resultRgccak_Tau0")
    prodYA1=resultRgccak_Tau0$Y[,1]%*%t(resultRgccak_Tau0$a[[1]])
    expect_true(all.equal(prodYA1,prodYA1_test))
})

test_that("test_rgcca_tau0ya2",{
    load(file="../results/resultRgccak_Tau0")
    prodYA2=resultRgccak_Tau0$Y[,2]%*%t(resultRgccak_Tau0$a[[2]])
    expect_true(all.equal(prodYA2,prodYA2_test))
})
test_that("test_rgcca_tau0ya3",{
    load(file="../results/resultRgccak_Tau0")
    prodYA3=resultRgccak_Tau0$Y[,3]%*%t(resultRgccak_Tau0$a[[3]])
    expect_true(all.equal(prodYA3,prodYA3_test))
})
# Tau optimal
#------------------
resultRgccak_TauOpt_test = rgccak(A, C, tau = rep("optimal",3), scheme = "factorial")

# test_that("test_rgccak_tauOpt",{
#       load(file="../results/resultRgccak_TauOpt")
#     # testing the same results than previously
#     expect_true(all.equal(resultRgccak_TauOpt_test$call$tau,resultRgccak_TauOpt[3:7]$tau))
# })

prodYA1_test=resultRgccak_TauOpt_test$Y[,1]%*%t(resultRgccak_TauOpt_test$a[[1]])
prodYA2_test=resultRgccak_TauOpt_test$Y[,2]%*%t(resultRgccak_TauOpt_test$a[[2]])
prodYA3_test=resultRgccak_TauOpt_test$Y[,3]%*%t(resultRgccak_TauOpt_test$a[[3]])

test_that("test_rgcca_tauOptya1",{
    load(file="../results/resultRgccak_TauOpt")
    prodYA1=resultRgccak_TauOpt$Y[,1]%*%t(resultRgccak_TauOpt$a[[1]])
    expect_true(all.equal(prodYA1,prodYA1_test))
})

test_that("test_rgcca_tauOptya2",{
    load(file="../results/resultRgccak_TauOpt")
    prodYA2=resultRgccak_TauOpt$Y[,2]%*%t(resultRgccak_TauOpt$a[[2]])
    expect_true(all.equal(prodYA2,prodYA2_test))
})
test_that("test_rgcca_tauOptya3",{
    load(file="../results/resultRgccak_TauOpt")
    prodYA3=resultRgccak_TauOpt$Y[,3]%*%t(resultRgccak_TauOpt$a[[3]])
    expect_true(all.equal(prodYA3,prodYA3_test))
})
# 
# set.seed(seed=2)
# A1=matrix(rnorm(500),10,50);rownames(A1)=paste0("S",1:10)
# set.seed(seed=3)
# A2=matrix(rnorm(300),10,30);rownames(A2)=paste0("S",1:10)
# set.seed(seed=4)
# A3=matrix(rnorm(20),10,2);rownames(A3)=paste0("S",1:10)
# C=matrix(1,3,3)-diag(1,3)
# blocks=list(A1=A1,A2=A2,A3=A3)
# resrgcca_k=RGCCA:::rgccak(A=blocks,C=C,scheme="factorial",tau=rep(1,3))
# #save(resrgcca,file="resRgcca_k_dual")
# load("../tests/results/resRgcca_k_dual")
# resrgcca_k$Y==resrgcca$Y
# 
# 

