# '# writeList test
# 
# '''
set.seed(42);X1=matrix(rnorm(35),7,5);X2=matrix(rnorm(28),7,4)
A=list(X1,X2)
writeList(A,wd=getwd())