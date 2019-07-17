# '#createNA3 test
# 
# '''
set.seed(42);X1=matrix(rnorm(35),7,5);X2=matrix(rnorm(28),7,4)
A=list(X1,X2)
# createNA
createNA(A,option="block",pNA=c(0.1,0.3))
createNA(A,option="rand",pNA=c(0.1))
