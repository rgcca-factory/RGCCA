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
library(RGCCA);
data(Russett);
namesFiles=dir("./R");
sapply(namesFiles,function(x){source(paste0("./R/",x))});

setwd("./Test/Results");
load(file="resultRgccak_Tau1");load(file="resultRgccak_Tau0");load(file="resultRgccak_TauOpt");
setwd("./../..");

# loading the new function output
 X_agric =as.matrix(Russett[,c("gini","farm","rent")]);
 X_ind = as.matrix(Russett[,c("gnpr","labo")]);
 X_polit = as.matrix(Russett[ , c("demostab", "dictator")]);
 A = list(X_agric, X_ind, X_polit);
 C = matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3);

T1=Sys.time();resultRgccak_Tau1_test = rgccak(A, C, tau = c(1, 1, 1), scheme = "factorial", scale = TRUE);T2=Sys.time();Tdiff_Tau1_test=T2-T1; #0.1548102 

T1=Sys.time();resultRgccak_Tau0_test= rgccak(A, C, tau = c(0, 0, 0), scheme = "factorial", scale = TRUE);T2=Sys.time();Tdiff_Tau0_test=T2-T1 ; # 0.008856535
T1=Sys.time();resultRgccak_TauOpt_test = rgccak(A, C, tau = rep("optimal",3), scheme = "factorial", scale = TRUE);T2=Sys.time();Tdiff_TauOpt_test=T2-T1; #0.008605003 secs

# testing the same results than previously
all.equal(resultRgccak_Tau1,resultRgccak_Tau1_test);
all.equal(resultRgccak_Tau0,resultRgccak_Tau0_test);
all.equal(resultRgccak_TauOpt,resultRgccak_TauOpt_test);

#---------------------------------------------------------------------
# Checking new functionalities : Comparison of new rgcca 
#---------------------------------------------------------------------

library(MASS)
X_agric =as.matrix(Russett[,c("gini","farm","rent")]);
X_ind = as.matrix(Russett[,c("gnpr","labo")]);
X_polit = as.matrix(Russett[ , c("demostab", "dictator")]);
X_ind[3,]=NA;
X_polit[1:2,]=NA;
A = list(X_agric, X_ind, X_polit);

C = matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3);
T1=Sys.time();resultRgccak_Tau1NA_test = rgccak(A, C, tau = c(1, 1, 1), scheme = "factorial", scale = TRUE);T2=Sys.time();Tdiff_Tau1=T2-T1 ;#0.002044916 secs/
T1=Sys.time();resultRgccak_Tau0NA_test = rgccak(A, C, tau = c(0, 0, 0), scheme = "factorial", scale = TRUE);T2=Sys.time();Tdiff_Tau0=T2-T1;# 0.001989126
T1=Sys.time();resultRgccak_TauOptNA_test =rgccak(A, C, tau = rep("optimal",3), scheme = "factorial", scale = TRUE);T2=Sys.time();Tdiff_TauOpt=T2-T1;#0.002758265

# ajout du parametre estimateNA=TRUE
A = list(X_agric, X_ind, X_polit);
Aref=lapply(A,scale)

A[[1]][1,1]=NA
A[[1]][3,1]=NA
A[[1]][2,2]=NA
A[[2]][4,1]=NA
A[[2]][5,1]=NA
A[[2]][6,2]=NA
C=matrix(1,3,3)-diag(3)
A=lapply(A,scale2)
rgccak(A,C,tau=c(1,1,1), scheme = "factorial", scale = TRUE,verbose=TRUE,na.rm=TRUE,tol=1e-7)
rgccak(A,C,tau=c(1,1,1), scheme = "factorial", scale = TRUE,verbose=TRUE,na.rm=TRUE,estimateNA="no")
rgccak(A,C,tau=c(1,1,1), scheme = "factorial", scale = TRUE,verbose=TRUE,na.rm=TRUE,estimateNA="first")
rgccak(A,C,tau=c(1,1,1), scheme = "factorial", scale = TRUE,verbose=TRUE,na.rm=TRUE,estimateNA="iterative")
rgccak(A,C,tau=c(1,1,1), scheme = "factorial", scale = TRUE,verbose=TRUE,na.rm=TRUE,estimateNA="superblock")

# Tests pour l'imputation dans RGCCAk

# Liste des RGCCA de reference:
C=matrix(1,3,3)-diag(3)
A = list(X_agric, X_ind, X_polit);
A1=lapply(A,scale2)
A2=A1
A2[[1]]=A1[[1]]/sqrt(3)
A2[[2]]=A1[[2]]/sqrt(2)
A2[[3]]=A1[[3]]/sqrt(2)


#----------------------
# Pour la fonction horst
#--------------------------
# on recupere la rgcca reference
refRgcca1=rgccak(A1,C,tau=c(1,1,1), scheme = "horst", verbose=TRUE,na.rm=TRUE,sameBlockWeight=FALSE,tol=1e-9)
refRgcca2=rgccak(A2,C,tau=c(1,1,1), scheme = "horst", verbose=TRUE,na.rm=TRUE,sameBlockWeight=TRUE,tol=1e-9)

# on cree les jeux de données avec valeur manquantes
A_miss0=createNA(A=Aref,option="ponc",pNA=0.1,nAllRespondants=10,output="list")$dat
    #==========================================================================
    # /!\ Remarque : Si on ne scale pas à la base le A, cela ne marche pas ! ! 
    #==========================================================================

# Si on n'a pas le same block weight :
A_miss1=lapply(A_miss0,scale2)
# on reproduit le traitement effectué dans RGCCAk si on a le scaling et le same block weight
A_miss2=A_miss1
A_miss2[[1]]=A_miss1[[1]]/sqrt(3)
A_miss2[[2]]=A_miss1[[2]]/sqrt(2)
A_miss2[[3]]=A_miss1[[3]]/sqrt(2)


res1=rgccak(A_miss1,C,tau=c(1,1,1), scheme = "horst", verbose=TRUE,na.rm=TRUE,estimateNA="lebrusquet",sameBlockWeight=FALSE,tol=1e-9)

nip1=rgccak(A_miss1,C,tau=c(1,1,1), scheme = "horst", verbose=TRUE,na.rm=TRUE,tol=1e-9)
nip2=rgccak(A_miss2,C,tau=c(1,1,1), scheme = "horst", verbose=TRUE,na.rm=TRUE,tol=1e-9)
res2=rgccak(A_miss2,C,tau=c(1,1,1), scheme = "horst", verbose=TRUE,na.rm=TRUE,estimateNA="lebrusquet",sameBlockWeight=TRUE,tol=1e-9)

plot(res1$crit,type="b", main="RGCCA criterion",col="blue")
lines(refRgcca1$crit,col="green")
lines(nip1$crit,col="red")

plot(res2$crit,type="b", main="RGCCA criterion",col="blue")
lines(refRgcca2$crit,col="green")
lines(nip2$crit,col="red")

plot(res2$crit)
lines(nip2$crit)
head(res1$A[[1]] )
head( A_miss1[[1]])
plot(A1[[1]][,1],res1$A[[1]][,1],col=1+is.na(A_miss0[[1]][,1]))
cor(res2$a[[1]],refRgcca2$a[[1]])
cor(nip2$a[[1]],refRgcca2$a[[1]])

cor(res2$Y[,1],refRgcca2$Y[,1])
cor(nip2$Y[,1],refRgcca2$Y[,1])

cor(res1$a[[1]],refRgcca$a[[1]])
cor(nip1$a[[1]],refRgcca$a[[1]])

#----------------------
# Pour une autre fonction 
#--------------------------
scheme="factorial"
# on recupere la rgcca reference
refRgcca1=rgccak(A1,C,tau=c(1,1,1), scheme = scheme, verbose=TRUE,na.rm=TRUE,sameBlockWeight=FALSE,tol=1e-9)
refRgcca2=rgccak(A2,C,tau=c(1,1,1), scheme = scheme, verbose=TRUE,na.rm=TRUE,sameBlockWeight=TRUE,tol=1e-9)

# on cree les jeux de données avec valeur manquantes
A_miss0=createNA(A=Aref,option="ponc",pNA=0.1,nAllRespondants=10,output="list")$dat

# Si on n'a pas le same block weight :
A_miss1=lapply(A_miss0,scale2)
# on reproduit le traitement effectué dans RGCCAk si on a le scaling et le same block weight
A_miss2=A_miss1
A_miss2[[1]]=A_miss1[[1]]/sqrt(3)
A_miss2[[2]]=A_miss1[[2]]/sqrt(2)
A_miss2[[3]]=A_miss1[[3]]/sqrt(2)

# si on ne scale pas à la base le A, cela ne marche pas ! ! 
res1=rgccak(A_miss1,C,tau=c(1,1,1), scheme = scheme, verbose=TRUE,na.rm=TRUE,estimateNA="lebrusquet",sameBlockWeight=FALSE,tol=1e-5)
nip1=rgccak(A_miss1,C,tau=c(1,1,1), scheme = scheme, verbose=TRUE,na.rm=TRUE,tol=1e-9)
nip2=rgccak(A_miss2,C,tau=c(1,1,1), scheme = scheme, verbose=TRUE,na.rm=TRUE,tol=1e-9)
res2=rgccak(A_miss2,C,tau=c(1,1,1), scheme = scheme, verbose=TRUE,na.rm=TRUE,estimateNA="lebrusquet",sameBlockWeight=TRUE,tol=1e-5)

plot(res1$crit,type="b", main="RGCCA criterion",col="blue")
lines(refRgcca1$crit,col="green")
lines(nip1$crit,col="red")

plot(res2$crit,type="b", main="RGCCA criterion",col="blue")
lines(refRgcca2$crit,col="green")
lines(nip2$crit,col="red")

plot(res2$crit)
lines(nip2$crit)
head(res1$A[[1]] )
head( A_miss1[[1]])
plot(A1[[1]][,1],res1$A[[1]][,1],col=1+is.na(A_miss0[[1]][,1]))
plot(A1[[1]][,2],res1$A[[1]][,2],col=1+is.na(A_miss0[[1]][,2]))

cor(res2$a[[1]],refRgcca2$a[[1]])
cor(nip2$a[[1]],refRgcca2$a[[1]])

cor(res2$Y[,1],refRgcca2$Y[,1])
cor(nip2$Y[,1],refRgcca2$Y[,1])

cor(res1$a[[1]],refRgcca$a[[1]])
cor(nip1$a[[1]],refRgcca$a[[1]])


