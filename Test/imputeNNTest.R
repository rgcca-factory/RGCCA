# '# Test imputeNN
# 
# '''
set.seed(42);X1=matrix(rnorm(50),10,5);
set.seed(22);X2=matrix(rnorm(40),10,4);
set.seed(2);X3=matrix(rnorm(70),10,7);
# usual test 
X1[1,]=NA
X2[7,1]=NA
X2[5,1]=NA
X3[3,]=NA
X3[4,]=NA
A=list(X1,X2,X3)
imputeNN(A,k=1,output="mean", klim=NULL,scale=TRUE,sameBlockWeight=TRUE)
  