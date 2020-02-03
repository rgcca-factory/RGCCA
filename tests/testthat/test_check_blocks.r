#'# imputeEM test

#'''
A=list(agri=Russett[,1:3],ind=Russett[,4:5],polit=Russett[,7:8])
Amiss=createNA(A,option="block",pNA=rep(0.1,3))$dat
C=matrix(1,3,3)- diag(3);tau=rep(1,3)
resEM=imputeEM(Amiss,C=C,tau=tau,ncomp=rep(2,3))
