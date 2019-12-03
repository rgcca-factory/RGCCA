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

test_that("test_rgcca_tau1",{
    load(file="../results/resultRgccak_Tau1")
    expect_true(all.equal(resultRgccak_Tau1_test,resultRgccak_Tau1))
})

load(file="../results/resultRgccak_Tau0")
test_that("test_rgccak_tau0",{
    resultRgccak_Tau0_test = rgccak(A, C, tau = c(0, 0, 0), scheme = "factorial")
     # testing the same results than previously
    expect_true(all.equal(resultRgccak_Tau0_test,resultRgccak_Tau0))
})

test_that("test_rgccak_tauOpt",{
      load(file="../results/resultRgccak_TauOpt")
    resultRgccak_TauOpt_test = rgccak(A, C, tau = rep("optimal",3), scheme = "factorial")
    # testing the same results than previously
    expect_true(all.equal(resultRgccak_TauOpt_test,resultRgccak_TauOpt))
})


