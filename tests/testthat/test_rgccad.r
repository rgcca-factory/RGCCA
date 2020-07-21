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

T1=Sys.time();resultRgcca_Tau1_test = rgccad(A, C, ncomp=rep(2,3),tau = c(1, 1, 1), scheme = "factorial", scale = TRUE,verbose=FALSE);T2=Sys.time();Tdiff_Tau1_test=T2-T1; #0.03056502
T1=Sys.time();resultRgcca_Tau0_test= rgccad(A, C, ncomp=rep(2,3),tau = c(0, 0, 0), scheme = "factorial", scale = TRUE,verbose=FALSE);T2=Sys.time();Tdiff_Tau0_test=T2-T1;  # 0.03978729
T1=Sys.time();resultRgcca_TauOpt_test = rgccad(A, C, ncomp=rep(2,3),tau = rep("optimal",3), scheme = "factorial", scale = TRUE,verbose=FALSE);T2=Sys.time();Tdiff_TauOpt_test=T2-T1 ;# 0.04306483   secs

names(A)=c("agri","ind","demo")
res = rgccad(A, C, ncomp=rep(2,3),tau = c(1, 0, 1), scheme = "factorial", scale = TRUE,verbose=FALSE)

# testing pca
res=rgccad(list(a1=A[[1]],b1=A[[1]]),C=matrix(c(0,1,1,0),2,2),tau=rep(1,2),ncomp=c(2,2))
test_that("test_rgcca_ave1",{
expect_true(round(res$AVE$AVE_X[[1]][1],digits=4)==0.7433)
})
test_that("test_rgcca_ave1",{
    expect_true(round(res$AVE$AVE_X[[1]][2],digits=4)==0.2371)
})

# testing the same results than previously
rgccaTest1=round(abs(cor(resultRgcca_Tau1$a[[1]][,1],resultRgcca_Tau1_test$a[[1]][,1])),digits=7)==1

test_that("test_rgcca_tau1",{
      expect_true(rgccaTest1)
})

# tests with different number of components
X_agric =as.matrix(Russett[,c("gini","farm","rent")]);
X_ind = as.matrix(Russett[,c("gnpr","labo")]);
X_polit = as.matrix(Russett[ , c("demostab","dictator")]);
A = list(X_agric, X_ind, X_polit);
res=rgccad(A, C, ncomp=c(2,2,1),tau = c(1, 1, 1), scheme = "factorial", scale = TRUE,verbose=FALSE)


res=rgccad(A, C, ncomp=c(2,2,2),tau = rep("optimal",3), scheme = "factorial",verbose=FALSE)


 #############
 # Example 1 #
 #############
 data(Russett)
 X_agric =as.matrix(Russett[,c("gini","farm","rent")])
 X_ind = as.matrix(Russett[,c("gnpr","labo")])
 X_polit = as.matrix(Russett[ , c("demostab", "dictator")])
 A = list(X_agric, X_ind, X_polit)
 #Define the design matrix (output = C) 
 C = matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3)
 result.rgcca = rgccad(A, C, tau = c(1, 1, 1), scheme = "factorial", scale = TRUE,verbose=FALSE)
 lab = as.vector(apply(Russett[, 9:11], 1, which.max))
 plot(result.rgcca$Y[[1]], result.rgcca$Y[[2]], col = "white", 
      xlab = "Y1 (Agric. inequality)", ylab = "Y2 (Industrial Development)")
 text(result.rgcca$Y[[1]], result.rgcca$Y[[2]], Russett[, 1], col = lab, cex = .7)

 ############################################
 # Example 2: RGCCA and multiple components #
 ############################################
 ############################
 # plot(y1, y2) for (RGCCA) #
 ############################
 result.rgcca = rgccad(A, C, tau = rep(1, 3), ncomp = c(2, 2, 1),
                      scheme = "factorial", verbose = FALSE)
layout(t(1:2))
 plot(result.rgcca$Y[[1]][, 1], result.rgcca$Y[[2]][, 1], col = "white", xlab = "Y1 (GE)", 
 ylab = "Y2 (CGH)", main = "Factorial plan of RGCCA")
 text(result.rgcca$Y[[1]][, 1], result.rgcca$Y[[2]][, 1], Russett[, 1], col = lab, cex = .6)
 plot(result.rgcca$Y[[1]][, 1], result.rgcca$Y[[1]][, 2], col = "white", xlab = "Y1 (GE)", 
     ylab = "Y2 (GE)", main = "Factorial plan of RGCCA")
 text(result.rgcca$Y[[1]][, 1], result.rgcca$Y[[1]][, 2], Russett[, 1], col = lab, cex = .6)


