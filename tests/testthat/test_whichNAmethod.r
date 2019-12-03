library(parallel)
library(FactoMineR)

# test on random data
#---------------------
set.seed(42);X1=matrix(rnorm(350),70,5);X2=matrix(rnorm(280),70,4);colnames(X1)=paste("A",1:5);colnames(X2)=paste("B",1:4);
rownames(X1)=rownames(X2)=paste("S",1:70)
A=list(X1,X2);


resultComparison=whichNAmethod(A,listMethods=c("nipals","mean"),patternNA=rep(0.1,2),seed=1:20)
resultComparison2=whichNAmethod(A,listMethods=c("nipals","mean"),patternNA=rep(0.1,2),seed=1:20)
test_that("whichNAmethod_1",{expect_true(all.equal(resultComparison,resultComparison2))})
