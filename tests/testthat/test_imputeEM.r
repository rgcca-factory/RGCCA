data(Russett)
A=list(agri=Russett[,1:3],ind=Russett[,4:5],polit=Russett[,7:8])
Amiss=createNA(A,option="ponc",pNA=rep(0.2,3))$dat
C=matrix(1,3,3)- diag(3);tau=rep(1,3)
resEM=imputeEM(A=Amiss,C=C,tau=tau,ncomp=rep(2,3))

plot(resEM$crit)
par(las=1)
par(mfrow=c(1,2))
plot(resEM$crit,pch=16,main="RGCCA criterion evolution",ylab="",xlab="Iteration")
plot(resEM$stab,pch=16,main="Stability",ylab="",xlab="Iteration")
