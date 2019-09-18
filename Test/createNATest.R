# '#createNA3 test
# 
# '''
set.seed(42);X1=matrix(rnorm(500),100,5);X2=matrix(rnorm(400),100,4)
rownames(X1)=rownames(X2)=paste("S",1:dim(X1)[1],sep="")
A=list(X1,X2)
# createNA
jdd=createNA(A,option="block",pNA=c(0.1,0.3))
sum(is.na(jdd$dat[[1]][,1]))
sum(is.na(jdd$dat[[2]][,1]))

jdd$dat[[1]][jdd$subjectKept,]
jdd$dat[[2]][jdd$subjectKept,]


res=createNA(A,option="block",pNA=c(0.5,0.5))

jdd=createNA(A,option="ponc",pNA=c(0.1,0.3))
jdd$dat[[1]][jdd$subjectKept,]
jdd$dat[[2]][jdd$subjectKept,]