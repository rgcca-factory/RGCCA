data(Russett)
A=list(agri=Russett[,1:3],ind=Russett[,4:5],polit=Russett[,7:8])
Amiss=A
Amiss[[1]][1:2,1:2]=NA
Amiss[[1]][3:4,2:3]=NA
Amiss[[1]][5,3]=NA
Amiss[[2]][1,1]=NA
Amiss[[2]][6,2]=NA
Amiss[[3]][7:8,]=NA
Amiss[[3]][9,2]=NA

C=matrix(1,3,3)- diag(3);tau=rep(1,3)
resEM=impute(blocks=Amiss,connection=C,tau=tau,ncomp=rep(1,3))


block=1
indmiss=is.na(Amiss[[block]])
plot(A[[block]][indmiss],resEM$A[[block]][indmiss],xlab="True values",ylab="Estimated values")
test_that("Checking RMSE value on an example",{expect_true(
round(sqrt(sum((A[[block]][indmiss]-resEM$A[[block]][indmiss])^2)/length(indmiss)),digits=6)==1.895868
)})

abline(a=0,b=1)

plot(resEM$crit)
par(las=1)
par(mfrow=c(1,2))
plot(resEM$crit,pch=16,main="RGCCA criterion evolution",ylab="",xlab="Iteration")
plot(resEM$stab,pch=16,main="Stability",ylab="",xlab="Iteration")

