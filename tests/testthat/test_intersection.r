# '# Test intersection
# 
# '''
set.seed(42);X1=matrix(rnorm(35),7,5);
set.seed(22);X2=matrix(rnorm(28),7,4);
set.seed(2);X3=matrix(rnorm(49),7,7);
# usual test 
X1[1,]=NA
X2[7,1]=NA
X2[5,1]=NA
A=list(X1,X2)
Ainter=intersection(A=A)
test_that("intersection_1",{expect_true(dim(Ainter[[1]])[1]==4)})
# too many subjects with missing values
X3[3,1:2]=NA
Ainter2=intersection(A=list(X1,X2,X3))
