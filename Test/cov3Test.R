'# cov3 test

'''
set.seed(42);X1=matrix(rnorm(35),7,5);
set.seed(22);X2=matrix(rnorm(35),7,5);

# Check old functionalities
cov3(X1,bias=FALSE)==cov(X1)
cov3(X1,X2,bias=FALSE)==cov(X1,X2)


# with missing values
X2[2:4,3:4]=NA
cov3(X2,bias=FALSE)==cov(X2)

# when one subject is missing
X1[1,]=NA
cov(X1)
cov3(X1,bias=FALSE)