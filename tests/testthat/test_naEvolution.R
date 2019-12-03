
set.seed(42);X1=matrix(rnorm(350),70,5);X2=matrix(rnorm(280),70,4)
A=list(X1,X2)
listResults=NULL
listResults=naEvolution(A=A,prctNA=c(0.1,0.2,0.3,0.4),listMethods=c("mean","complete"))
test_that("naEvolution_1",{expect_true(!is.null(listResults))})