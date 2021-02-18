# '# Test intersection_list
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
Ainter=intersection_list(A=A)

X1[1,]=NA
X2[7,1]=NA
X2[5,1]=NA
A2=lapply(A,scale)
Ainter=intersection_list(A=A2)

test_that("intersection_list_1",{expect_true(dim(Ainter[[1]])[1]==4)})
# too many subjects with missing values
X3[3,1:2]=NA
expect_error(
  intersection_list(A=list(X1,X2,X3)), 
  "Less than 3 subjects have no missing values, choose another
               missing value handling method or work on your dataset.", 
  fixed=TRUE
)
