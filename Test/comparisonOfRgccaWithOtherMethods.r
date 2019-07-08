# setwd("/home/caroline.peltier/Bureau/RGCCA")
rm(list=ls())
library(RGCCA)
library(MASS)
library(nipals)
namesFiles=dir("./R")
# loading functions in R directory
sapply(namesFiles,function(x){source(paste0("./R/",x))})

#-----------------------------------------------
# Loading example
#-----------------------------------------------
data(Russett)
X_agric =as.matrix(Russett[,c("gini","farm","rent")])
X_ind = as.matrix(Russett[,c("gnpr","labo")])
X_polit = as.matrix(Russett[ , c("demostab", "dictator")])
refData= list(X_agric, X_ind, X_polit)

#-----------------------------------------------
# find back unscaled PCA without superblock
#-----------------------------------------------
C=matrix(1,2,2);diag(C)=0
A=list(refData[[1]],refData[[1]])
resRgcca=rgcca(A,C,scale=FALSE,ncomp=rep(2,2))
resPca=prcomp((A2[[1]]),scale=FALSE)
cor(resPca$x[,1:2],resRgcca$Y[[1]]) # composante 2 inversée
cor(resRgcca$a[[1]][,1],resPca$rotation[,1])
cor(resRgcca$a[[1]][,2],resPca$rotation[,2])

#-----------------------------------------------
# Find back unscaled PCA without superblock
#-----------------------------------------------
C=matrix(1,2,2);diag(C)=0
A=list(refData[[1]],refData[[1]])
resRgcca=rgcca(A,C,scale=TRUE,ncomp=rep(2,2))
resPca=prcomp((A2[[1]]),scale=TRUE)
cor(resPca$x[,1:2],resRgcca$Y[[1]]) # composante 2 inversée
cor(resRgcca$a[[1]][,1],resPca$rotation[,1])
cor(resRgcca$a[[1]][,2],resPca$rotation[,2])

#-----------------------------------------------
# Finding back the results of RGCCA with SGCCA
#-----------------------------------------------
A=list(refData[[1]],refData[[2]])
resRgcca=rgcca(A,C,scale=TRUE,ncomp=rep(2,2))
resSgcca=sgcca(A,C,scale=TRUE,ncomp=rep(2,2),c1=rep(1,2))
cor(resSgcca$Y[[1]],resRgcca$Y[[1]])


resSgcca=sgcca(A2,C,scale=TRUE,ncomp=rep(2,2),c1=rep(1,2),na.rm=FALSE)
resSgcca$a[[1]]
resSgcca=sgcca(A2,C,scale=TRUE,ncomp=rep(2,2),c1=rep(1/sqrt(2),2),na.rm=FALSE)
resSgcca$a[[1]]
resSgcca=sgcca(A2,C,scale=TRUE,ncomp=rep(2,2),c1=rep(1/sqrt(3),2),na.rm=FALSE)
resSgcca$a[[1]]
resSgcca=sgcca(A2,C,scale=TRUE,ncomp=rep(2,2),c1=rep(1/sqrt(4),2),na.rm=FALSE)
resSgcca$a[[1]]
resSgcca=sgcca(A2,C,scale=TRUE,ncomp=rep(2,2),c1=rep(1/sqrt(5),2),na.rm=FALSE)
resSgcca$a[[1]]
#--------------------
# Finding back PCA with superblock
#-------------------------

X = Russett[, 1:8]
fit.pca =prcomp(X, scale = TRUE)

L = c(as.list(X), list(X))
L = lapply(L, as.matrix)
names(L) = NULL

C = matrix(0, 9, 9)
C[, 9] = C[, 9] = 1
C[9, ] = C[9, ] = 1
diag(C) = 0

fit.rgcca = rgcca(A=L, tau = c(rep(0, 8), 1),
                    ncomp = c(rep(1, 8), 2),
                    scheme = "factorial",
                    scale = TRUE, init = "random",
                    bias = FALSE, tol = 1e-16)
y = fit.rgcca$Y[[9]]
a = fit.rgcca$a[[9]]
cor(a[,1:2],fit.pca$rotation[,1:2])
beta= apply(scale2(X), 2, function(x) lm(x~y)$coefficients[2])
cor(y,fit.pca$x[,1:2])
#---------------------
# retrouver les valeurs manquantes
#---------------------
X = Russett[, 1:8]
fit.pca =prcomp(X, scale = TRUE)

L = c(as.list(X), list(X))
L = lapply(L, as.matrix)
names(L) = NULL

C = matrix(0, 9, 9)
C[, 9] = C[, 9] = 1
diag(C) = 0

fit.rgcca = rgcca(L, tau = c(rep(0, 8), 1),
                  ncomp = c(rep(1, 8), 2),
                  scheme = "factorial",
                  scale = TRUE, init = "random",
                  bias = FALSE, tol = 1e-16)

y = fit.rgcca$Y[[9]]
a = apply(scale2(X), 2, function(x) lm(x~y)$coefficients[2])

# Données manquantes 
indNA = cbind(sample(1:NROW(X), 20, replace = TRUE), sample(1:NCOL(X), 20, replace = TRUE))
indNA = indNA[!duplicated(indNA), ]
XNA = X
XNA[indNA] = NA

X1NA = XNA
m = apply(X1NA, 2, function(x) mean(x, na.rm = TRUE))
X1NA[indNA] = m[indNA[, 2]]

for (i in 1:10){
  L = c(as.list(as.data.frame(X1NA)), list(X1NA))
  fit.rgcca = rgcca(L, tau = c(rep(0, 8), 0),
                    ncomp = c(rep(1, 8), 1),
                    scheme = "factorial",
                    scale = TRUE, init = "random",
                    verbose = FALSE, tol = 1e-16)
  
  y = fit.rgcca$Y[[9]]
  a = apply(scale2(X1NA, bias = TRUE), 2, function(x) lm(x~y)$coefficients[2])
  
  X2NA = scale2(X1NA, bias = TRUE)
  moy = matrix(attr(X2NA, "scaled:center"),
               nrow = NROW(X2NA), 8, byrow = TRUE)
  stdev = matrix(attr(X2NA, "scaled:scale"),
                 nrow = NROW(X2NA), 8, byrow = TRUE)
  Xhat = (y%*%t(a))*stdev + moy
  X1NA[indNA] = Xhat[indNA]

}
#- em pca
library()
fit.PCA.EM = imputePCA(XNA, ncp = 1, method = "EM",maxiter=11,threshold=1e-16)
# my implementation
L = as.list(as.data.frame(XNA)); C=matrix(1,length(L),length(L));diag(C)=0
res=imputeSB(L,C=C,tau=rep(0,length(L)),niter=9)$A
do.call(cbind,res)[indNA]
fit.PCA.EM $completeObs[indNA]
#--------------------
# retrouver l'AFM avec un superblock
#-------------------------
X = Russett[, 1:8]
library(FactoMineR)
fit.MFA= MFA(X,group=c(3,5))

L = c(list(X[,1:3]),list(X[,4:5]),list(X[,6:8]), list(X))
names(L) = NULL

C = matrix(0, 4, 4)
C[, 4] = C[, 4] = 1
C[4, ] = C[4, ] = 1
diag(C) = 0

fit.rgcca = rgcca1(L, tau = c(rep(0, 3), 0),
                  ncomp = c(rep(1, 3), 1),C=C,
                  scheme = "factorial",
                  scale = TRUE, init = "svd",
                  verbose = FALSE, tol = 1e-16,sameBlockWeight=T)
y = fit.rgcca$Y[[4]]
a = fit.rgcca$a[[4]]

cor(fit.MFA$ind$coord[,1],fit.rgcca[["Y"]][[4]][,1])

fit.MFA$global.pca$var$coord[,1:2]
fit.MFA$global.pca$ind$coord[,1:2][1:5,]

y[1:5,]
#--------------------
# retrouver la pls ??
#-----------------------
A3=list(refData[[1]],refData[[2]][,1])
A4=lapply(A3,scale,scale=FALSE)
resplsr=plsr(A4[[2]]~A4[[1]],validation="none")
resRgcca2=rgcca(A4,C,scale=FALSE,ncomp=c(2,1))
resplsr$Yscores
resplsr$Yloadings
resRgcca2$Y

#-------------------------------------
# retrouver l'analyse de correlation ??
#--------------------------------------
rescancor=cancor(t(refData[[1]]),t(refData[[2]]))
resRgcca1=rgcca(list(refData[[1]],refData[[2]]),C=C,tau=rep(0,2),ncomp=rep(5,2))


#---------------------------------------------------------
# Find the binding of PCA components with the superbloc
#---------------------------------------------------------
refData = list(X_agric, X_ind, X_polit)
nBloc=length(refData)
superbloc=do.call(cbind,refData)
C=matrix(0,nBloc+1,nBloc+1);C[,nBloc+1]=1;C[nBloc+1,]=1;C[nBloc+1,nBloc+1]=0
A=c(refData,list(superbloc))
refDataS=lapply(refData[1:nBloc],scale)
nBloc=length(refDataS)
superblocS=do.call(cbind,refDataS)
As=c(refDataS,list(superblocS))

# Concatenation 
### scale = TRUE
res=rgcca(As,C=C,tau=c(rep(1,nBloc),1),scheme="factorial",tol=1e-16,scale=TRUE,sameBlockWeight = FALSE)
cor( res$a[[nBloc+1]][1:3] ,res$a[[1]])   # on retrouve exactement les mêmes axes
cor( res$a[[nBloc+1]][4:5] ,res$a[[2]]) 
### scale = FALSE
res=rgcca(A,C=C,tau=rep(1,nBloc+1),scheme="factorial",tol=1e-16,sameBlockWeight = FALSE, scale = FALSE)
fit.pca = prcomp(do.call(cbind, res$Y[1:3]), scale = FALSE) # getting the supercomponent
cor(fit.pca$x[, 1], res$Y[[4]]) #Y of superblock is correlated to the supercomponent of the three Y by block
### option superblock
res=rgcca(A,C=C,tau=c(rep(1,nBloc),1),scheme="factorial",tol=1e-16,scale=TRUE,sameBlockWeight = TRUE,superblock=TRUE)
cor( res$a[[nBloc+1]][1:3] ,res$a[[1]])   # cor=1 : concatenation des vecteurs


# Finding back the supercomponent of Y
### scale = TRUE -> on ne retrouve pas ! !
res=rgcca(A,C=C,tau=rep(1,nBloc+1),scheme="factorial",tol=1e-16,sameBlockWeight = FALSE, scale = TRUE)
respca=prcomp(Reduce(cbind,res$Y[1:3]),scale=TRUE)
cor(respca$x[,1],res$Y[[1]])

### scale = FALSE : on retrouve
res=rgcca(A,C=C,tau=rep(1,nBloc+1),scheme="factorial",tol=1e-16,sameBlockWeight = FALSE, scale = FALSE)
fit.pca = prcomp(do.call(cbind, res$Y[1:3]), scale = FALSE) # PCA of first block: "supercomponent"
cor(fit.pca$x[, 1], res$Y[[4]])

resS=rgcca(As,C=C,tau=rep(1,nBloc+1),scheme="factorial",tol=1e-16,sameBlockWeight = FALSE, scale = FALSE)
fit.pcaS = prcomp(do.call(cbind, resS$Y[1:3]), scale = FALSE) # PCA of first block: "supercomponent"
cor(fit.pcaS$x[, 1], resS$Y[[4]])

### option superblock
#### scale = FALSE
res=rgcca(A,C=C,tau=c(rep(1,nBloc),1),scheme="factorial",tol=1e-16,scale=FALSE,sameBlockWeight = TRUE,superblock=TRUE)
fit.pcaS = prcomp(do.call(cbind, resS$Y[1:3]), scale = FALSE) 
cor(fit.pcaS$x[, 1], resS$Y[[4]])
#### scale = TRUE
res=rgcca(A,C=C,tau=c(rep(1,nBloc),1),scheme="factorial",tol=1e-16,scale=TRUE,sameBlockWeight = TRUE,superblock=TRUE)
fit.pcaS = prcomp(do.call(cbind, resS$Y[1:3]), scale = TRUE) 
cor(fit.pcaS$x[, 1], resS$Y[[4]])

 resPcaSuperblock=prcomp(Reduce(cbind,refData[1:3]),scale=TRUE)
 cor(resPcaSuperblock$x[,1],respca$x[,1])
 
 # Verifier que si on mélange les variables par bloc on obtient la même composante pour le superbloc
 res=rgcca(A,C=C,tau=c(rep(1,nBloc),1),scheme="factorial",tol=1e-16,scale=TRUE,sameBlockWeight = FALSE)
 refDataS2=c(list(cbind(refDataS[[1]],refDataS[[2]])),(refDataS[3]))
 nBloc2=length(refDataS2)
 superbloc2=do.call(cbind,refDataS2)
 C2=matrix(0,nBloc2+1,nBloc2+1);C2[,nBloc2+1]=1;C2[nBloc2+1,]=1;C2[nBloc2+1,nBloc2+1]=0
 A2=c(refDataS2,list(superbloc2))
 res2=rgcca(A2,C=C2,tau=rep(1,nBloc2+1),scheme="factorial",tol=1e-16,sameBlockWeight = FALSE)
 cor(res$Y[[nBloc+1]], res2$Y[[nBloc2+1]]) # On retrouve exactement la même chose
 cor(res$a[[nBloc+1]], res2$a[[nBloc2+1]]) # ici aussi
 
 
 # selection de variables
 refData=refData[1:3]
superblock=do.call(cbind,refData)
superblockS=apply(superblock,2,scale) #necessite que toutes les variables aient une variance non nulle
refDataS=lapply(refData,scale)
C=matrix(0,length(refDataS)+1,length(refDataS)+1);C[length(refDataS)+1,]=1;C[,(length(refDataS)+1)]=1;diag(C)=0
resRgcca1=rgccaNa(c(refDataS,list(superblockS)),C=C,tau=rep(1,length(refDataS)+1),na.impute="none",ncomp=rep(2,length(refDataS)+1),sameBlockWeight = F,scheme="factorial")
mat=cor(superblockS)
matDist=1-abs(mat)
resHclust=hclust(dist(matDist) )
plot(resHclust)
plot(rev(resHclust$height))
nGroup=8
varGroup=cutree(hclust(dist(matDist)),nGroup)
listRes=list()
listRes=lapply(1:nGroup,function(i){listRes[[i]]=superblockS[,varGroup==i]})
C2=matrix(0,length(listRes)+1,length(listRes)+1);C2[length(listRes)+1,]=1;C2[,(length(listRes)+1)]=1;diag(C2)=0
listResF=c(listRes,list(superblockS));lapply(listResF,dim)
resRgcca2=rgccaNa(listResF,C=C2,tau=rep(1,length(listRes)+1),na.impute="none",ncomp=rep(2,length(listRes)+1),sameBlockWeight=FALSE,scheme="factorial")
cor(resRgcca2$a[[length(listRes)+1]][,1],resRgcca1$a[[length(refDataS)+1]][,1]) # validé !
#----------Gestion des valeurs manquantes dans cette demarche
setwd("/home/caroline.peltier/Bureau/EtudeNA/Datasets/bloc10/1")
source("/home/caroline.peltier/Bureau/NewFunctionsRGCCA/scale3.r")

testData=readDataset()[1:3]
superblock=do.call(cbind,testData[1:3])
superblockS=scale3(superblock)
refDataS=scale3(do.call(cbind,refData[1:3])) # on scale toute la liste de reference
mat=cor(superblockS,use = "pairwise") # on calcule les correlations
matDist=1-abs(mat)
resHclust=hclust(dist(matDist) )
plot(resHclust)
plot(rev(resHclust$height))
nGroup=4
varGroup=cutree(hclust(dist(matDist),),nGroup) # on obtient une classification de variables
listRes=list()
listRes=lapply(1:nGroup,function(i){listRes[[i]]=superblockS[,varGroup==i]}) # on decoupe suivant cette classification
lapply(listRes,dim)
# la classification permet elle d'avoir des blocs sans "lignes manquantes" ?
lapply(listRes,function(x){
          res=rep(TRUE,dim(x)[1])
          for(i in 1:dim(x)[1])
          {
              if(sum(is.na(x[i,]))==dim(x)[2]){res[i]=FALSE}
          }
          return(res)
        })
C=matrix(1,length(listRes)+1,length(listRes)+1);C[length(listRes)+1,]=1;C[,length(listRes)+1]=1;diag(C)=0
resRgcca=rgccaNa(c(listRes,list(superblockS)),C=C,tau=rep(1,length(listRes)+1),na.impute="iterative",ncomp=rep(2,length(listRes)+1),sameBlockWeight=FALSE,scheme="factorial")
Anew=(Reduce(cbind,resRgcca$A))
W=is.na(do.call(cbind,testData))
Anew[W]
refDataS[W]

  