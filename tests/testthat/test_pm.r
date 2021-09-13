# '# Test pm
#
# '''
test_that("test_pm",{
    set.seed(42); A=matrix(rnorm(15),3,5)
    set.seed(34); B=matrix(rnorm(20),5,4)
    # loading the new function output
    pmToGet=A%*%B;
    pmRes=pm(A,B);
    # testing the same results than previously
    expect_true(all.equal(pmToGet,pmRes))
})

test_that("test_pmNA",{
    # loading the new function output
    set.seed(42); A=matrix(rnorm(15),3,5)
    set.seed(34); B=matrix(rnorm(20),5,4)
    A[3,5]=0
    B[1,1]=0
    T1=Sys.time();pmToGet=A%*%B;T2=Sys.time(); Tdiff1=T2-T1
    A[3,5]=NA
    B[1,1]=NA
    T1=Sys.time();pmRes=pm(A,B);T2=Sys.time();Tdiff2=T2-T1
    expect_true(all.equal(pmToGet,pmRes))
})
