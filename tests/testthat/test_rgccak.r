# '# Test rgccak
#
# '''
# setwd("./tests/testthat")
library(MASS)
data(Russett)
X_agric =as.matrix(Russett[, c("gini","farm","rent")]);
X_ind = as.matrix(Russett[, c("gnpr","labo")]);
X_polit = as.matrix(Russett[ , c("demostab", "dictator")]);
A = list(X_agric, X_ind, X_polit);
C = matrix(c(0, 0, 1,
             0, 0, 1,
             1, 1, 0), 3, 3);

resultRgccak_Tau1_test = rgccak(A, C, tau = c(1, 1, 1), scheme = "factorial")
prodYA1_test = resultRgccak_Tau1_test$Y[, 1]%*%t(resultRgccak_Tau1_test$a[[1]])
prodYA2_test = resultRgccak_Tau1_test$Y[, 2]%*%t(resultRgccak_Tau1_test$a[[2]])
prodYA3_test = resultRgccak_Tau1_test$Y[, 3]%*%t(resultRgccak_Tau1_test$a[[3]])

# test_that("test_rgcca_tau1",{
#     load(file="../results/resultRgccak_Tau1")
#     expect_true(all.equal(resultRgccak_Tau1_test[3:7],resultRgccak_Tau1[3:7]))
# })

test_that("test_rgcca_tau1ya1",{
    load(file = "../results/resultRgccak_Tau1")
    prodYA1 = resultRgccak_Tau1$Y[, 1]%*%t(resultRgccak_Tau1$a[[1]])
     expect_true(all.equal(abs(prodYA1), abs(prodYA1_test)))
  })

test_that("test_rgcca_tau1ya2",{
    load(file = "../results/resultRgccak_Tau1")
    prodYA2 = resultRgccak_Tau1$Y[, 2]%*%t(resultRgccak_Tau1$a[[2]])
    expect_true(all.equal(prodYA2, prodYA2_test))
})

test_that("test_rgcca_tau1ya3",{
    load(file = "../results/resultRgccak_Tau1")
    prodYA3 = resultRgccak_Tau1$Y[, 3]%*%t(resultRgccak_Tau1$a[[3]])
    expect_true(all.equal(abs(prodYA3), abs(prodYA3_test)))
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
prodYA1_test = resultRgccak_Tau0_test$Y[, 1]%*%t(resultRgccak_Tau0_test$a[[1]])
prodYA2_test = resultRgccak_Tau0_test$Y[, 2]%*%t(resultRgccak_Tau0_test$a[[2]])
prodYA3_test = resultRgccak_Tau0_test$Y[, 3]%*%t(resultRgccak_Tau0_test$a[[3]])

test_that("test_rgcca_tau0ya1",{
    load(file = "../results/resultRgccak_Tau0")
    prodYA1 = resultRgccak_Tau0$Y[, 1]%*%t(resultRgccak_Tau0$a[[1]])
    expect_true(all.equal(abs(prodYA1), abs(prodYA1_test)))
})

test_that("test_rgcca_tau0ya2",{
    load(file = "../results/resultRgccak_Tau0")
    prodYA2 = resultRgccak_Tau0$Y[, 2]%*%t(resultRgccak_Tau0$a[[2]])
    expect_true(all.equal(abs(prodYA2), abs(prodYA2_test)))
})

test_that("test_rgcca_tau0ya3",{
    load(file = "../results/resultRgccak_Tau0")
    prodYA3=resultRgccak_Tau0$Y[, 3]%*%t(resultRgccak_Tau0$a[[3]])
    expect_true(all.equal(abs(prodYA3), abs(prodYA3_test)))
})

# Tau optimal
#------------------
resultRgccak_TauOpt_test = rgccak(A, C, tau = rep("optimal",3),
                                  scheme = "factorial")

prodYA1_test = resultRgccak_TauOpt_test$Y[, 1]%*%
  t(resultRgccak_TauOpt_test$a[[1]])
prodYA2_test = resultRgccak_TauOpt_test$Y[, 2]%*%
  t(resultRgccak_TauOpt_test$a[[2]])
prodYA3_test = resultRgccak_TauOpt_test$Y[, 3]%*%
  t(resultRgccak_TauOpt_test$a[[3]])

test_that("test_rgcca_tauOptya1",{
    load(file = "../results/resultRgccak_TauOpt")
    prodYA1 = resultRgccak_TauOpt$Y[, 1]%*%t(resultRgccak_TauOpt$a[[1]])
    expect_true(all.equal(abs(prodYA1), abs(prodYA1_test)))
})

test_that("test_rgcca_tauOptya2",{
    load(file = "../results/resultRgccak_TauOpt")
    prodYA2 = resultRgccak_TauOpt$Y[, 2]%*%t(resultRgccak_TauOpt$a[[2]])
    expect_true(all.equal(abs(prodYA2), abs(prodYA2_test)))
})

test_that("test_rgcca_tauOptya3",{
    load(file = "../results/resultRgccak_TauOpt")
    prodYA3=resultRgccak_TauOpt$Y[, 3]%*%t(resultRgccak_TauOpt$a[[3]])
    expect_true(all.equal(abs(prodYA3), abs(prodYA3_test)))
})


# Cas primal with missing values
#----------------------------------------
X_agric =as.matrix(Russett[,c("gini", "farm", "rent")]);
X_ind = as.matrix(Russett[,c("gnpr", "labo")]);
X_polit = as.matrix(Russett[ , c("demostab", "dictator")]);
A = list(X_agric, X_ind, X_polit);
A[[3]][1:3, 1:2]=NA
A[[3]][5, 1]=NA
res_tau1_primal = rgccak(A, C = 1-diag(3),scheme="factorial")
res_tau0_primal = rgccak(A, C = 1-diag(3),scheme="factorial",
                         tau = c(0.5, 0.5, 0.5))

# Cas dual
#------------
data(Russett)
X_agric =as.matrix(Russett[1:9, 1])
X_ind = as.matrix(Russett[1:9, 2:11])
A=list(ag = X_agric, ind = X_ind)
A[[2]][1:3, 1:3]=NA
j = 2
t(!is.na(A[[j]]))%*%(!is.na(A[[j]]))
res_tau0_dual=rgccak(A, C = 1 - diag(2), tau=c(0.5, 0.5))
res_tau1_dual=rgccak(A, C = 1 - diag(2), tau=c(1, 1))

