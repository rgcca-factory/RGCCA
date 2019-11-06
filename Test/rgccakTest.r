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
rgccak(A,C,tau=c(1,1,1), scheme = "factorial", verbose=TRUE,na.rm=TRUE,tol=1e-7)
rgccak(A,C,tau=c(1,1,1), scheme = "factorial", verbose=TRUE,na.rm=TRUE,estimateNA="no")
rgccak(A,C,tau=c(1,1,1), scheme = "factorial", verbose=TRUE,na.rm=TRUE,estimateNA="first")
rgccak(A,C,tau=c(1,1,1), scheme = "factorial", verbose=TRUE,na.rm=TRUE,estimateNA="iterative")
rgccak(A,C,tau=c(1,1,1), scheme = "factorial", verbose=TRUE,na.rm=TRUE,estimateNA="superblock")


#==========================================================================
# Tests pour la fonction leb imputeInRgcca
#==========================================================================

#----------------------
# Pour Russetts
#----------------------
X_agric =as.matrix(Russett[,c("gini","farm","rent")]);
X_ind = as.matrix(Russett[,c("gnpr","labo")]);
X_polit = as.matrix(Russett[ , c("demostab", "dictator")]);
C=matrix(1,3,3)-diag(3)
A = list(X_agric, X_ind, X_polit);
A1=lapply(A,scale2)
A2=A1
A2[[1]]=A1[[1]]/sqrt(3)
A2[[2]]=A1[[2]]/sqrt(2)
A2[[3]]=A1[[3]]/sqrt(2)
# on recupere la rgcca reference
refRgcca1=rgccak(A1,C,tau=c(1,1,1), scheme = "horst",init="random", verbose=TRUE,na.rm=TRUE,sameBlockWeight=FALSE,tol=1e-14)
refRgcca2=rgccak(A2,C,tau=c(1,1,1), scheme = "horst", verbose=TRUE,na.rm=TRUE,sameBlockWeight=FALSE,tol=1e-9)

#perm=sgcca.permute.crit(A1,c1s=matrix(1,nrow=2,ncol=3),nperm=1000,scheme="horst",plot=TRUE)
nTests=20
correlR2=correlR1=correlN2=correlN1=pb1=pb2=rep(NA,nTests)
corLL2=corNipals2=nLL2=nNipals2=corLL=corNipals=nLL=nNipals=rep(NA,nTests)
for(i in 1:nTests)
{
    # creation des jeux de données sur le jdd initial
    A_miss0=createNA(A=A,option="ponc",pNA=0.005,nAllRespondants=10,output="list")$dat
    # centrage reduction du jeu de données obtenu
    A_miss1=lapply(A_miss0,scale2)
    # on reproduit le traitement effectué dans RGCCAk si on a le scaling et le same block weight
    A_miss2=A_miss1
    # division par rapport à j
    A_miss2[[1]]=A_miss1[[1]]/sqrt(3)
    A_miss2[[2]]=A_miss1[[2]]/sqrt(2)
    A_miss2[[3]]=A_miss1[[3]]/sqrt(2)
    res1=rgccak(A_miss1,C,tau=c(1,1,1), scheme = "horst", verbose=TRUE,na.rm=TRUE,estimateNA="lebrusquet",sameBlockWeight=FALSE,tol=1e-5)
    nip1=rgccak(A_miss1,C,tau=c(1,1,1), scheme = "horst", verbose=FALSE,na.rm=TRUE,tol=1e-5)
    nip2=rgccak(A_miss2,C,tau=c(1,1,1), scheme = "horst", verbose=FALSE,na.rm=TRUE,tol=1e-5)
    res2=rgccak(A_miss2,C,tau=c(1,1,1), scheme = "horst", verbose=TRUE,na.rm=TRUE,estimateNA="lebrusquet",sameBlockWeight=TRUE,tol=1e-5)
    
    pb1[i]=sum(diff(res1$crit)<0)
    pb2[i]=sum(diff(res2$crit)<0)
    
    #plot(res1$crit,type="b", main="RGCCA criterion",col="blue")
    #lines(refRgcca1$crit,col="green")
    #lines(nip1$crit,col="red")
    
    #plot(res2$crit,type="b", main="RGCCA criterion",col="blue")
    #lines(refRgcca2$crit,col="green")
    #lines(nip2$crit,col="red")
    
    #head(res1$A[[1]] )
    #head( A_miss1[[1]])
   # plot(A1[[1]][,1],res1$A[[1]][,1],col=1+is.na(A_miss0[[1]][,1]))
    abline(a=0,b=1)
    plot(A1[[1]][,2],res1$A[[1]][,2],col=1+is.na(A_miss0[[1]][,2]))
    abline(a=0,b=1)
    
    mean(res1$A[[1]][,1])
    sd(res1$A[[1]][,1])*sqrt(46/47)
    
    correlR2[i]=cor(res2$a[[1]],refRgcca2$a[[1]])
    correlN2[i]=cor(nip2$a[[1]],refRgcca2$a[[1]])
    correlR1[i]=cor(res1$a[[1]],refRgcca1$a[[1]])
    correlN1[i]=cor(nip1$a[[1]],refRgcca1$a[[1]])
    
    nLL[i]= min(sum((res1$a[[1]]-refRgcca1$a[[1]])^2),sum((res1$a[[1]]+refRgcca1$a[[1]])^2))
    nNipals[i]= min(sum((nip1$a[[1]]-refRgcca1$a[[1]])^2),sum((nip1$a[[1]]+refRgcca1$a[[1]])^2))
    nLL2[i]= min(sum((res2$a[[1]]-refRgcca2$a[[1]])^2),sum((res2$a[[1]]+refRgcca2$a[[1]])^2))
    nNipals2[i]= min(sum((nip2$a[[1]]-refRgcca2$a[[1]])^2),sum((nip2$a[[1]]+refRgcca2$a[[1]])^2))
    corLL2[i]=abs(cor(    res2$a[[1]],    refRgcca2$a[[1]]))
    corNipals2[i]=abs(cor(    nip2$a[[1]],    refRgcca2$a[[1]]))
}

ds=data.frame(c(rep("Nipals",nTests),rep("Leb",nTests)),c(log(nNipals2),log(nLL2)))
colnames(ds)=c("meth","rep")
library(ggplot2)
p <- ggplot(ds, aes(meth, rep))
p + geom_boxplot()

summary(corLL)
summary(corNipals)
summary(nLL)
summary(nNipals)
ds=data.frame(c(rep("Nipals",nTests),rep("Leb",nTests)),c(log(nNipals),log(nLL)))
colnames(ds)=c("meth","rep")
library(ggplot2)
p <- ggplot(ds, aes(meth, rep))
p + geom_boxplot()

#--------------------------------
# Pour un jeu de données simulé
#--------------------------------
nTests=50
n=100
Z=rnorm(n)
sdNoise=0.1
X1=cbind(Z+rnorm(n=n,sd=sdNoise),-0.5*(Z)+3+rnorm(n,sd=sdNoise),Z+rnorm(n,sdNoise))
X2=cbind(2-rnorm(m=rnorm(n),sd=sdNoise,n)*Z,3*Z+rnorm(n=n))
A=list(X1,X2) # jeu de données initial
nBlock=length(A)
A1=lapply(A,scale2,bias=TRUE) # jeu de données centré réduit
A2=lapply(A1,function(x){return(x/sqrt(ncol(x)))}) #divisé par racine de p

sigma = matrix(c(1,0.2,0.8,
                 0.2,1,0.7,
                 0.8,0.7,1),3,3, byrow = TRUE)

X = mvrnorm(50, rep(0,3),Sigma = sigma)
v1 <- matrix(c(rep(.5,25),rep(0,75)),ncol=1)
v2 <- matrix(c(rep(1,25),rep(0,25)),ncol=1)
v3 <- matrix(c(rep(.5,25),rep(0,175)),ncol=1)
c = 0.5
x1 <- X[,1]%*%t(v1) + c*matrix(rnorm(50*100),ncol=100)
x2 <- X[,2]%*%t(v2) + c*matrix(rnorm(50*50),ncol=50)
x3 <- X[,3]%*%t(v3) + c*matrix(rnorm(50*200),ncol=200)
xlist <- lapply(list(x1, x2, x3),scale2)

#A1=xlist;nBlock=length(A1)

refRgcca1=rgccak(A=A1,C=matrix(1,nBlock,nBlock)-diag(nBlock),tau=rep(1,nBlock),sameBlockWeight=FALSE,init="svd",scheme="horst",tol=1e-5)
#refRgcca2=rgccak(A=A2,C=matrix(1,nBlock,nBlock)-diag(nBlock),tau=rep(1,nBlock),sameBlockWeight=TRUE,init="svd",scheme="horst",tol=1e-5)
corLL2=corNipals2=nLL2=nNipals2=corLL=corNipals=nLL=nNipals=rep(NA,nTests)
for(i in 1:nTests)
{
    A=list(X1,X2)
    #A=list(x1,x2,x3)
    A_miss=createNA(A=A,option="ponc",pNA=0.005,nAllRespondants=10,output="list")$dat
    A1_miss=lapply(A_miss,scale2,bias=TRUE)
    A2_miss=lapply(A1_miss,function(x){return(x/sqrt(NCOL(x)))})
    nip1=rgccak(A=A1_miss,C=matrix(1,nBlock,nBlock)-diag(nBlock),tau=rep(1,nBlock),sameBlockWeight=FALSE,init="svd",scheme="horst",verbose=FALSE,tol=1e-10)
    res1=rgccak(A=A1_miss,C=matrix(1,nBlock,nBlock)-diag(nBlock),tau=rep(1,nBlock),sameBlockWeight=FALSE,init="svd",estimateNA="lebrusquet",scheme="horst",verbose=TRUE,tol=1e-10,na.rm=TRUE)
    
  #  nip2=rgccak(A=A2_miss,C=matrix(1,nBlock,nBlock)-diag(nBlock),tau=rep(1,nBlock),sameBlockWeight=TRUE,init="svd",scheme="horst",verbose=FALSE,tol=1e-5)
  #  res2=rgccak(A=A2_miss,C=matrix(1,nBlock,nBlock)-diag(nBlock),tau=rep(1,nBlock),sameBlockWeight=TRUE,init="svd",estimateNA="lebrusquet",scheme="horst",verbose=TRUE,tol=1e-5,na.rm=TRUE)
       
    corLL[i]=abs(cor(    res1$a[[1]],    refRgcca1$a[[1]]))
    corNipals[i]=abs(cor(    nip1$a[[1]],    refRgcca1$a[[1]]))
    nLL[i]= min(sum((res1$a[[1]]-refRgcca1$a[[1]])^2),sum((res1$a[[1]]+refRgcca1$a[[1]])^2))
    nNipals[i]= min(sum((nip1$a[[1]]-refRgcca1$a[[1]])^2),sum((nip1$a[[1]]+refRgcca1$a[[1]])^2))
    
    # corLL2[i]=abs(cor(    res2$a[[1]],    refRgcca2$a[[1]]))
    # corNipals2[i]=abs(cor(    nip2$a[[1]],    refRgcca2$a[[1]]))
    # nLL2[i]= min(sum((res1$a[[1]]-refRgcca2$a[[1]])^2),sum((res2$a[[1]]+refRgcca2$a[[1]])^2))
    # nNipals2[i]= min(sum((nip2$a[[1]]-refRgcca2$a[[1]])^2),sum((nip2$a[[1]]+refRgcca2$a[[1]])^2))
}

summary(corLL)
summary(corNipals)
summary(nLL)
summary(nNipals)
ds=data.frame(c(rep("Nipals",nTests),rep("Leb",nTests)),c(log(nNipals),log(nLL)))
colnames(ds)=c("meth","rep")
library(ggplot2)
p <- ggplot(ds, aes(meth, rep))
p + geom_boxplot()







