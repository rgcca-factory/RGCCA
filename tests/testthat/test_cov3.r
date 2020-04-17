#'# cov3 test

#'''
set.seed(42);X1=matrix(rnorm(35),7,5);
set.seed(22);X2=matrix(rnorm(35),7,5);

# Check old functionalities
test_that("cov3_1",{expect_true(sum(cov3(X1,bias=FALSE)==cov(X1))==25)})
test_that("cov3_2",{expect_true(sum(cov3(X1,X2,bias=FALSE)==cov(X1,X2))==25)})


# with missing values
X2[2:4,3:4]=NA
test_that("cov3_3",{expect_true(sum(cov3(X2,bias=FALSE)==cov(X2),na.rm=TRUE)==9)})

# when one subject is missing
X1[1,]=NA
cov(X1)
cov3(X1,bias=FALSE)

