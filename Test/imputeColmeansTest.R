'# Test imputeColmeans

'''
A=matrix(rnorm(15),3,5);A[3,5]=NA
B=matrix(rnorm(15),3,5);A[3,]=NA
imputeColmeans(list(A,B))

