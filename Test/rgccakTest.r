# '# Test RGCCAk
# 
# '''

#setwd("/home/caroline.peltier/Bureau/RGCCA")

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
# T1=Sys.time();resultRgccak_Tau1 = RGCCA::rgccak(A, C, tau = c(1, 1, 1), scheme = "factorial", scale = TRUE);T2=Sys.time();Tdiff_Tau1=T2-T1 #0.002044916 secs/
# T1=Sys.time();resultRgccak_Tau0 = RGCCA::rgccak(A, C, tau = c(0, 0, 0), scheme = "factorial", scale = TRUE);T2=Sys.time();Tdiff_Tau0=T2-T1# 0.001989126
# T1=Sys.time();resultRgccak_TauOpt = RGCCA::rgccak(A, C, tau = rep("optimal",3), scheme = "factorial", scale = TRUE);T2=Sys.time();Tdiff_TauOpt=T2-T1#0.002758265
# 
# save(file="resultRgccak_Tau1",resultRgccak_Tau1);save(file="resultRgccak_Tau0",resultRgccak_Tau0);save(file="resultRgccak_TauOpt",resultRgccak_TauOpt)
# setwd("./../..")


# loading the reference output
library(RGCCA)
data(Russett)
namesFiles=dir("./R")
sapply(namesFiles,function(x){source(paste0("./R/",x))})

setwd("./Test/Results")
load(file="resultRgccak_Tau1");load(file="resultRgccak_Tau0");load(file="resultRgccak_TauOpt")
setwd("./../..")

# loading the new function output


 X_agric =as.matrix(Russett[,c("gini","farm","rent")])
 X_ind = as.matrix(Russett[,c("gnpr","labo")])
 X_polit = as.matrix(Russett[ , c("demostab", "dictator")])
 A = list(X_agric, X_ind, X_polit)
 C = matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3)

T1=Sys.time();resultRgccak_Tau1_test = rgccak(A, C, tau = c(1, 1, 1), scheme = "factorial", scale = TRUE);T2=Sys.time();Tdiff_Tau1_test=T2-T1 #0.1548102 
T1=Sys.time();resultRgccak_Tau0_test= rgccak(A, C, tau = c(0, 0, 0), scheme = "factorial", scale = TRUE);T2=Sys.time();Tdiff_Tau0_test=T2-T1  # 0.008856535
T1=Sys.time();resultRgccak_TauOpt_test = rgccak(A, C, tau = rep("optimal",3), scheme = "factorial", scale = TRUE);T2=Sys.time();Tdiff_TauOpt_test=T2-T1 #0.008605003 secs

# testing the same results than previously
all.equal(resultRgccak_Tau1,resultRgccak_Tau1_test)
all.equal(resultRgccak_Tau0,resultRgccak_Tau0_test)
all.equal(resultRgccak_TauOpt,resultRgccak_TauOpt_test)

#---------------------------------------------------------------------
# Checking new functionalities : Comparison of new rgcca 
#---------------------------------------------------------------------


X_agric =as.matrix(Russett[,c("gini","farm","rent")])
X_ind = as.matrix(Russett[,c("gnpr","labo")])
X_ind[3,]=NA
X_polit = as.matrix(Russett[ , c("demostab", "dictator")])
X_polit[1:2,]=NA
A = list(X_agric, X_ind, X_polit)

C = matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3)
T1=Sys.time();resultRgccak_Tau1NA_test = rgccak(A, C, tau = c(1, 1, 1), scheme = "factorial", scale = TRUE);T2=Sys.time();Tdiff_Tau1=T2-T1 #0.002044916 secs/
T1=Sys.time();resultRgccak_Tau0NA_test = rgccak(A, C, tau = c(0, 0, 0), scheme = "factorial", scale = TRUE);T2=Sys.time();Tdiff_Tau0=T2-T1# 0.001989126
T1=Sys.time();resultRgccak_TauOptNA_test =rgccak(A, C, tau = rep("optimal",3), scheme = "factorial", scale = TRUE);T2=Sys.time();Tdiff_TauOpt=T2-T1#0.002758265



