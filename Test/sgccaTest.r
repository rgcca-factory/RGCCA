# '# Test SGCCA
# 
# '''



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
# C = matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3)
# T1=Sys.time();resultSgcca = RGCCA::sgcca(A, C,ncomp=rep(2,3), scheme = "factorial", scale = TRUE);T2=Sys.time();Tdiff_Tau1=T2-T1 #0.04458332 secs/
# 
# save(file="resultSgcca",resultSgcca)
# setwd("./../..")


# loading the reference output
setwd("./Test/Results")
load(file="resultSgcca");
setwd("./../..")

# loading the new function output
data(Russett)
X_agric =as.matrix(Russett[,c("gini","farm","rent")]);
X_ind = as.matrix(Russett[,c("gnpr","labo")]);
X_polit = as.matrix(Russett[ , c("demostab", "dictator")]);
A = list(X_agric, X_ind, X_polit);
C = matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3);

namesFiles=dir("./R");
sapply(namesFiles,function(x){source(paste0("./R/",x))});

T1=Sys.time();resultSgcca_test = sgcca(A, C, ncomp=rep(2,3), scheme = "factorial", scale = TRUE,verbose=FALSE);T2=Sys.time();Tdiff_Tau1_test=T2-T1; #0.03056502 

# testing the same results than previously
rgccaTest1=all.equal(resultSgcca,resultSgcca_test);
#---------------------------------------------------------------------
# Checking new functionalities : Comparison of new rgcca 
#---------------------------------------------------------------------
#  Missing values

data(Russett);
X_agric =as.matrix(Russett[,c("gini","farm","rent")]);
X_ind = as.matrix(Russett[,c("gnpr","labo")]);
X_ind[3,]=NA;
X_polit = as.matrix(Russett[ , c("demostab", "dictator")]);
X_polit[1:2,]=NA;
A = list(X_agric, X_ind, X_polit);
C = matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3);
T1=Sys.time();resultSgccaNA_test = sgcca(A, ncomp=rep(2,3),C,  scheme = "factorial", scale = TRUE,verbose=FALSE);T2=Sys.time();Tdiff_Tau1=T2-T1; #0.002044916 secs/

