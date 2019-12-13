#'# determine_patternNA

#'''
X1=matrix(rnorm(150),30,5)
X2=matrix(rnorm(150),30,5)
X3=matrix(rnorm(150),30,5)
X1[1:10,1]=NA
X1[2,2]=NA
X1[11,]=NA
X2[1:10,1]=NA
A=list(bloc1=X1,bloc2=X2,bloc3=X3)
res=determine_patternNA(A)
# Check old functionalities
test_that("determine_patternNA_1",{expect_true(sum(res$pctNA$bloc3)==0)})
test_that("determine_patternNA_2",{expect_true(length(res$pctNAbyBlock)==3)})
test_that("determine_patternNA_3",{expect_true(dim(res$completeSubjectByBlock)[2]==3)})
