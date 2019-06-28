'# Test RGCCAk

'''
library(MASS)
#library(RGCCA)
wd="/home/caroline.peltier/Bureau/RGCCAtoPull"
setwd(wd)



#---------------------------------------------------------------------
# Checking usual functionalities : Comparison of new rgcca with the packag for different tau
#---------------------------------------------------------------------
# # saving the reference output
# setwd("./Test/Results")
# data(Russett)
# X_agric =as.matrix(Russett[,c("gini","farm","rent")])
# X_ind = as.matrix(Russett[,c("gnpr","labo")])
# X_polit = as.matrix(Russett[ , c("demostab", "dictator")])
# A = list(X_agric, X_ind, X_polit)
# 
# C = matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3)
# T1=Sys.time();resultRgcca_Tau1 = RGCCA::rgcca(A, C,ncomp=rep(2,3), tau = c(1, 1, 1), scheme = "factorial", scale = TRUE);T2=Sys.time();Tdiff_Tau1=T2-T1 #0.04458332 secs/
# T1=Sys.time();resultRgcca_Tau0 = RGCCA::rgcca(A, C,ncomp=rep(2,3),tau = c(0, 0, 0), scheme = "factorial", scale = TRUE);T2=Sys.time();Tdiff_Tau0=T2-T1# 0.03887248
# T1=Sys.time();resultRgcca_TauOpt = RGCCA::rgcca(A, C,ncomp=rep(2,3), tau = rep("optimal",3), scheme = "factorial", scale = TRUE);T2=Sys.time();Tdiff_TauOpt=T2-T1#0.03995538
# 
# save(file="resultRgcca_Tau1",resultRgcca_Tau1);save(file="resultRgcca_Tau0",resultRgcca_Tau0);save(file="resultRgcca_TauOpt",resultRgcca_TauOpt)
# setwd("./../..")


# loading the reference output
setwd("./Test/Results")
load(file="resultRgcca_Tau1");load(file="resultRgcca_Tau0");load(file="resultRgcca_TauOpt")
setwd("./../..")

# loading the new function output
data(Russett)
X_agric =as.matrix(Russett[,c("gini","farm","rent")])
X_ind = as.matrix(Russett[,c("gnpr","labo")])
X_polit = as.matrix(Russett[ , c("demostab", "dictator")])
A = list(X_agric, X_ind, X_polit)
C = matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3)

source("rgcca.r")
source("rgccak.r")
source("pm.r")
source("defl.select.r")
source("initsvd.r")

T1=Sys.time();resultRgcca_Tau1_test = rgcca(A, C, ncomp=rep(2,3),tau = c(1, 1, 1), scheme = "factorial", scale = TRUE);T2=Sys.time();Tdiff_Tau1_test=T2-T1 #0.03056502 
T1=Sys.time();resultRgcca_Tau0_test= rgcca(A, C, ncomp=rep(2,3),tau = c(0, 0, 0), scheme = "factorial", scale = TRUE);T2=Sys.time();Tdiff_Tau0_test=T2-T1  # 0.03978729
T1=Sys.time();resultRgcca_TauOpt_test = rgcca(A, C, ncomp=rep(2,3),tau = rep("optimal",3), scheme = "factorial", scale = TRUE);T2=Sys.time();Tdiff_TauOpt_test=T2-T1 # 0.04306483   secs

# testing the same results than previously
all.equal(resultRgcca_Tau1,resultRgcca_Tau1_test)
all.equal(resultRgcca_Tau0,resultRgcca_Tau0_test)
all.equal(resultRgcca_TauOpt,resultRgcca_TauOpt_test) 
#---------------------------------------------------------------------
# Checking new functionalities : Comparison of new rgcca 
#---------------------------------------------------------------------
#  Missing values
setwd("./Test/Results")
data(Russett)
X_agric =as.matrix(Russett[,c("gini","farm","rent")])
X_ind = as.matrix(Russett[,c("gnpr","labo")])
X_ind[3,]=NA
X_polit = as.matrix(Russett[ , c("demostab", "dictator")])
X_polit[1:2,]=NA
A = list(X_agric, X_ind, X_polit)
library(nipals)
C = matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3)
T1=Sys.time();resultRgcca_Tau1NA_test = rgcca(A, ncomp=rep(2,3),C, tau = c(1, 1, 1), scheme = "factorial", scale = TRUE);T2=Sys.time();Tdiff_Tau1=T2-T1 #0.002044916 secs/
T1=Sys.time();resultRgcca_Tau0NA_test = rgcca(A, ncomp=rep(2,3),C, tau = c(0, 0, 0), scheme = "factorial", scale = TRUE);T2=Sys.time();Tdiff_Tau0=T2-T1# 0.001989126
T1=Sys.time();resultRgcca_TauOptNA_test =rgcca(A, C,ncomp=rep(2,3), tau = rep("optimal",3), scheme = "factorial", scale = TRUE);T2=Sys.time();Tdiff_TauOpt=T2-T1#0.002758265


# scale=FALSE gave false results: this test attests the difference between the results of previous and new RGCCA
data(Russett)
X_agric =as.matrix(Russett[,c("gini","farm","rent")])
X_ind = as.matrix(Russett[,c("gnpr","labo")])
X_polit = as.matrix(Russett[ , c("demostab", "dictator")])
A = list(X_agric, X_ind, X_polit)
C = matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3)
resultRgcca_nonScale_test = rgcca(A, C,ncomp=rep(2,3), tau = c(1, 1, 1), scheme = "factorial", scale = FALSE)
resultRgcca_nonScale  = RGCCA::rgcca(A, ncomp=rep(2,3),C, tau = c(1, 1, 1), scheme = "factorial", scale = FALSE)
all.equal(resultRgcca_nonScale ,resultRgcca_nonScale_test[1:10])


