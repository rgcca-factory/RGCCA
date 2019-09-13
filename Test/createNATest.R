# '#createNA3 test
# 
# '''
set.seed(42);X1=matrix(rnorm(350),70,5);X2=matrix(rnorm(280),70,4)
rownames(X1)=rownames(X2)=paste("S",1:dim(X1)[1],sep="")
A=list(X1,X2)
# createNA
jdd=createNA(A,option="block",pNA=c(0.1,0.3))

createNA(A,option="rand",pNA=c(0.1))

res=createNA(A,option="block",pNA=c(0.5,0.5))

