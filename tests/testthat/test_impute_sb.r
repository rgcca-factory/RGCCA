data(Russett)
A=list(agri=Russett[,1:5],ind=Russett[,4:5],polit=Russett[,7:8])
# MFA et rgcca 
# library(FactoMineR)
# nvar = sapply(A, NCOL)
# resmfa=MFA(base=do.call(cbind,A), group=nvar,  type=rep("s",length(nvar)))
# fit.rgcca=rgcca(A,scheme="factorial",ncomp=6,scale=TRUE,scale_block="lambda1",superblock=TRUE)
# resmfa$ind$coord[1:5,1:3]
# fit.rgcca$Y[[4]][1:5,1:3] #ok

#Creation of a dataset with NA
Amiss=A
Amiss[[1]][1:2,1:2]=NA
Amiss[[1]][3:4,2:3]=NA
Amiss[[1]][5,3]=NA
Amiss[[2]][1,1]=NA
Amiss[[2]][6,2]=NA
Amiss[[3]][7:8,]=NA
Amiss[[3]][9,2]=NA
resEM_sb=impute_sb(blocks =Amiss,scheme="factorial",scale=TRUE,scale_block="lambda1",ncomp=rep(1,3),bias=FALSE)

resEM_sb2=impute_sb(blocks=Amiss,scheme="factorial",scale=TRUE,scale_block=TRUE,ncomp=rep(1,3),bias=FALSE)
# TODO
#resEM_sb3=impute_sb(A=Amiss,scheme="factorial",scale=FALSE,scale_block="lambda1",ncomp=rep(1,3),bias=FALSE)

# Comparaison with MFA (ni =5 and naxis= 1)
library(missMDA)
nvar = sapply(Amiss, NCOL)
A_mfa=imputeMFA(X=do.call(cbind,Amiss), group=nvar, ncp = 1, type=rep("s",length(nvar)), method = "em",maxiter=5,threshold=1e-2)
resEM_sb=impute_sb(blocks=Amiss,scheme="factorial",scale=TRUE,scale_block="lambda1",ncomp=rep(1,3),bias=FALSE,ni=5)

test_that("Checking that imputeMFA is found back with rgcca (and 1 component)",{expect_true(
  sum(abs(round(A_mfa$completeObs-Reduce(cbind,resEM_sb$A),digits=5)))==0
)})

A_mfa2=imputeMFA(X=do.call(cbind,Amiss), group=nvar,ncp = 2, type=rep("s",length(nvar)), method = "em",maxiter=5,threshold=1e-2)
resEM_sb2=impute_sb(blocks=Amiss,scheme="factorial",scale=TRUE,scale_block="lambda1",naxis=2,bias=FALSE,ni=5)

test_that("Checking that imputeMFA is found back with rgcca (and 2 components)",{expect_true(
  sum(abs(round(A_mfa2$completeObs-Reduce(cbind,resEM_sb2$A),digits=5)))==0
)})

# NCOMP = 2

# TODO
#resEM_sb=impute_sb(A=Amiss,scheme="factorial",scale=TRUE,scale_block="lambda1",ncomp=2,naxis=2,bias=FALSE,ni=5)
