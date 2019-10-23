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
# X_polit = as.matrix(Russett[ , c("demostab", "dictatur")])
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
X_agric =as.matrix(Russett[,c("gini","farm","rent")]);
X_ind = as.matrix(Russett[,c("gnpr","labo")]);
X_polit = as.matrix(Russett[ , c("demostab", "dictator")]);
A = list(X_agric, X_ind, X_polit);
C = matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3);

T1=Sys.time();resultRgcca_Tau1_test = rgcca(A, C, ncomp=rep(2,3),tau = c(1, 1, 1), scheme = "factorial", scale = TRUE,verbose=FALSE);T2=Sys.time();Tdiff_Tau1_test=T2-T1; #0.03056502 
T1=Sys.time();resultRgcca_Tau0_test= rgcca(A, C, ncomp=rep(2,3),tau = c(0, 0, 0), scheme = "factorial", scale = TRUE,verbose=FALSE);T2=Sys.time();Tdiff_Tau0_test=T2-T1;  # 0.03978729
T1=Sys.time();resultRgcca_TauOpt_test = rgcca(A, C, ncomp=rep(2,3),tau = rep("optimal",3), scheme = "factorial", scale = TRUE,verbose=FALSE);T2=Sys.time();Tdiff_TauOpt_test=T2-T1 ;# 0.04306483   secs

# testing the same results than previously
rgccaTest1=all.equal(resultRgcca_Tau1,resultRgcca_Tau1_test);
rgccaTest2=all.equal(resultRgcca_Tau0,resultRgcca_Tau0_test);
rgccaTest3=all.equal(resultRgcca_TauOpt,resultRgcca_TauOpt_test) ;
rgccaTest1
rgccaTest2
rgccaTest3
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
T1=Sys.time();resultRgcca_Tau1NA_test = rgcca(A, ncomp=rep(2,3),C, tau = c(1, 1, 1), scheme = "factorial", scale = TRUE,verbose=FALSE);T2=Sys.time();Tdiff_Tau1=T2-T1; #0.002044916 secs/
T1=Sys.time();resultRgcca_Tau0NA_test = rgcca(A, ncomp=rep(2,3),C, tau = c(0, 0, 0), scheme = "factorial", scale = TRUE,verbose=FALSE);T2=Sys.time();Tdiff_Tau0=T2-T1;# 0.001989126

#T1=Sys.time();resultRgcca_TauOptNA_test =rgcca(A, C,ncomp=rep(2,3), tau = rep("optimal",3), scheme = "factorial", scale = TRUE);T2=Sys.time();Tdiff_TauOpt=T2-T1#0.002758265


# scale=FALSE gave false results: this test attests the difference between the results of previous and new RGCCA
data(Russett);
X_agric =as.matrix(Russett[,c("gini","farm","rent")]);
X_ind = as.matrix(Russett[,c("gnpr","labo")]);
X_polit = as.matrix(Russett[ , c("demostab", "dictatur")]);
A = list(X_agric, X_ind, X_polit);
C = matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3);
resultRgcca_nonScale_test = rgcca(A, C,ncomp=rep(2,3), tau = c(1, 1, 1), scheme = "factorial", scale = FALSE,verbose=FALSE);
resultRgcca_nonScale  = RGCCA::rgcca(A, ncomp=rep(2,3),C, tau = c(1, 1, 1), scheme = "factorial", scale = FALSE,verbose=FALSE);
rgccaTest4=all.equal(resultRgcca_nonScale ,resultRgcca_nonScale_test[1:10]);

# other test
data(Russett);
X_agric =as.matrix(Russett[,c("gini","farm","rent")]);
X_ind = as.matrix(Russett[,c("gnpr","labo")]);
X_polit = as.matrix(Russett[ , c("demostab", "dictator")]);
X_agric[c(2,4),]=NA;
X_ind[1,]=NA;
X_polit[5,1]=NA;
A = list(X_agric, X_ind, X_polit);
C = matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3);
rgcca(A=A,ncomp=rep(2,3),scale=TRUE,sameBlockWeight=TRUE,tau= c(1, 1, 1),verbose=FALSE);


X_agric =as.matrix(Russett[,c("gini","farm","rent")])
X_ind = as.matrix(Russett[,c("gnpr","labo")])
X_polit = as.matrix(Russett[ , c("demostab", "dictator")])
A = list(X_agric, X_ind, X_polit)
C = matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3)
result.rgcca = RGCCA::rgcca(A, C, tau = c(1, 1, 1), scheme = "factorial", scale = TRUE)
head(as.matrix(scale(X_agric)*sqrt(47/46)/sqrt(3))%*%result.rgcca$a[[1]])
head(result.rgcca$Y[[1]])
# on retrouve bien le resultat...
# avec RGCCA maintenant...
namesFiles=dir("./R");
sapply(namesFiles,function(x){source(paste0("./R/",x))});
result.rgcca.test = rgcca(A, C, tau = c(1, 1, 1), scheme = "factorial", scale = TRUE,returnA=TRUE)
head(as.matrix(scale(X_agric)*sqrt(47/46)/sqrt(3))%*%result.rgcca.test$a[[1]])
head(result.rgcca.test$Y[[1]])

# On teste le parametre estimateNA=TRUE pour les options suivantes : 
# -scale=TRUE, tau=1, dual=FALSE et ncomp=rep(1)
A = list(X_agric, X_ind, X_polit);
Aref=lapply(A,scale)
A[[1]][1,2]=NA
A[[2]][c(3,4,5),]=NA
A[[3]][c(6:10),]=NA

C=matrix(1,3,3)-diag(3)
result.rgcca.test1 = rgcca(A, C, tau = c(1, 1, 1), ncomp=rep(1,3),scheme = "factorial", scale = TRUE,returnA=TRUE,estimateNA="iterative",sameBlockWeight=FALSE,tol=1e-20)
result.rgcca.test2 = rgcca(A, C=, tau = c(1, 1, 1), ncomp=rep(1,3),scheme = "factorial", scale = TRUE,returnA=TRUE,estimateNA="first",sameBlockWeight=FALSE,tol=1e-20)
result.rgcca.test3 = rgcca(A, C="superblock", tau = c(1, 1, 1), ncomp=rep(1,3),scheme = "factorial", scale = TRUE,returnA=TRUE,estimateNA="superblock",sameBlockWeight=FALSE,tol=1e-20)
result.rgcca.test1 = rgcca(A, C, tau = c(1, 1, 1), ncomp=rep(1,3),scheme = "factorial", scale = TRUE,returnA=TRUE,estimateNA="iterative",sameBlockWeight=FALSE,tol=1e-20,verbose=TRUE)
result.rgcca.test1 = rgcca(A, C, tau = c(1, 1, 1), ncomp=rep(1,3),scheme = "factorial", scale = TRUE,returnA=TRUE,estimateNA="iterative",sameBlockWeight=FALSE,tol=1e-6,verbose=TRUE)

plot(result.rgcca.test2$crit)
result.rgcca.test1$a[[1]]
result.rgcca.test2$a[[1]] # les valeurs manquantes explosent sans contrainte
result.rgcca.ref $a[[1]]
result.rgcca.ref = rgcca(Aref, C, tau = c(1, 1, 1), ncomp=rep(1,3),scheme = "factorial", scale = TRUE,returnA=TRUE,estimateNA="no",sameBlockWeight=FALSE,tol=1e-20)
result.rgcca.test$imputedA[[1]][1,2]
Aref[[2]][c(3,4,5),]
result.rgcca.test$imputedA[[2]][3:5,]
Aref[[3]][c(6:10),]
result.rgcca.test$imputedA[[3]][c(6:10),]

result.rgcca.nipals = rgcca(A, C, tau = c(1, 1, 1), ncomp=rep(1,3),scheme = "factorial", scale = TRUE,returnA=TRUE,estimateNA=FALSE,na.rm=TRUE)

data.frame(Reduce("rbind",result.rgcca.ref$a),
Reduce("rbind",result.rgcca.nipals$a),
Reduce("rbind",result.rgcca.test$a))
