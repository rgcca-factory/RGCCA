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
# T1=Sys.time();resultSgccak = RGCCA::sgccak(A, C, scheme = "factorial", scale = TRUE);T2=Sys.time();Tdiff_Tau1=T2-T1 #0.002044916 secs/
 
# save(file="resultSgccak",resultSgccak);
# setwd("./../..")


# loading the reference output
library(RGCCA);
data(Russett);
namesFiles=dir("./R");
sapply(namesFiles,function(x){source(paste0("./R/",x))});

setwd("./Test/Results");
load(file="resultSgccak")
setwd("./../..");

# loading the new function output


 X_agric =as.matrix(Russett[,c("gini","farm","rent")]);
 X_ind = as.matrix(Russett[,c("gnpr","labo")]);
 X_polit = as.matrix(Russett[ , c("demostab", "dictator")]);
 A = list(X_agric, X_ind, X_polit);
 C = matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3);

T1=Sys.time();resultSgccak_test = sgccak(A, C, scheme = "factorial", scale = TRUE);T2=Sys.time();Tdiff_Tau1_test=T2-T1; #0.1548102 

# testing the same results than previously
all.equal(resultSgccak,resultSgccak_test);

#---------------------------------------------------------------------
# Checking new functionalities : Comparison of new sgccak with missing values
#---------------------------------------------------------------------
X_agric =as.matrix(Russett[,c("gini","farm","rent")]);
X_ind = as.matrix(Russett[,c("gnpr","labo")]);
X_ind[3,]=NA;
X_polit = as.matrix(Russett[ , c("demostab", "dictator")]);
X_polit[1:2,]=NA;
A = list(X_agric, X_ind, X_polit);

C = matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3);
T1=Sys.time();resultSgccak_NA_test = sgccak(A, C,  scheme = "factorial", scale = TRUE);T2=Sys.time();Tdiff_Tau1=T2-T1 ;#0.002044916 secs/



