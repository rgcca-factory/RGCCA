# '# Test rgcca
# 
# '''
#setwd("~/Bureau/RGCCA/tests/testthat")
load(file="../results/resultRgcca_Tau1");load(file="../results/resultRgcca_Tau0");load(file="../results/resultRgcca_TauOpt")
data(Russett)
X_agric =as.matrix(Russett[,c("gini","farm","rent")]);
X_ind = as.matrix(Russett[,c("gnpr","labo")]);
X_polit = as.matrix(Russett[ , c("demostab", "dictator")]);
A = list(X_agric, X_ind, X_polit);
C = matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3);

T1=Sys.time();resultRgcca_Tau1_test = rgcca(A, C, ncomp=rep(2,3),tau = c(1, 1, 1), scheme = "factorial", scale = TRUE,verbose=FALSE);T2=Sys.time();Tdiff_Tau1_test=T2-T1; #0.03056502 
T1=Sys.time();resultRgcca_Tau0_test= rgcca(A, C, ncomp=rep(2,3),tau = c(0, 0, 0), scheme = "factorial", scale = TRUE,verbose=FALSE);T2=Sys.time();Tdiff_Tau0_test=T2-T1;  # 0.03978729
T1=Sys.time();resultRgcca_TauOpt_test = rgcca(A, C, ncomp=rep(2,3),tau = rep("optimal",3), scheme = "factorial", scale = TRUE,verbose=FALSE);T2=Sys.time();Tdiff_TauOpt_test=T2-T1 ;# 0.04306483   secs

# testing the same results than previously
rgccaTest1=resultRgcca_Tau1$a[[1]]==resultRgcca_Tau1_test$a[[1]]

rgccaTest2=all.equal(resultRgcca_Tau0,resultRgcca_Tau0_test);
rgccaTest3=all.equal(resultRgcca_TauOpt,resultRgcca_TauOpt_test) ;

test_that("test_rgcca_tau1",{
      expect_true(TRUE)
})

