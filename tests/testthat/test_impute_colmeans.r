# '# Test impute_colmeans
# 
# '''
A=matrix(rnorm(15),3,5);A[3,5]=NA
B=matrix(rnorm(15),3,5);B[3,]=NA
C=impute_colmeans(list(A,B))

test_that("impute_colmeans_1",{expect_true(C[[1]][3,5]==mean(A[,5],na.rm=TRUE))})
