library(parallel)


# test on random data
#---------------------
set.seed(42);X1=matrix(rnorm(350),70,5);X2=matrix(rnorm(280),70,4);colnames(X1)=paste("A",1:5);colnames(X2)=paste("B",1:4);
rownames(X1)=rownames(X2)=paste("S",1:70)
A=list(X1,X2);
names(A)=c("bloc1","bloc2")


resultComparison=whichNAmethod(A,listMethods=c("nipals","mean"),patternNA=rep(0.1,2),seed=1:20)
resultComparison2=whichNAmethod(A,listMethods=c("nipals","mean"),patternNA=rep(0.1,2),seed=1:20)
test_that("whichNAmethod_1",{expect_true(all.equal(resultComparison,resultComparison2))})


set.seed(42);X1=matrix(rnorm(350),5,70);X2=matrix(rnorm(350),5,70);colnames(X1)=paste("A",1:70);colnames(X2)=paste("B",1:70);
rownames(X1)=rownames(X2)=paste("S",1:5)
A=list(X1,X2);
names(A)=c("bloc1","bloc2")

resultComparison=whichNAmethod(A,ncomp=rep(2,4),listMethods=c("nipals","mean"),patternNA=rep(0.1,2),seed=1:20)
resultComparison2=whichNAmethod(A,listMethods=c("nipals","mean"),patternNA=rep(0.1,2),seed=1:20)
resultComparison3=whichNAmethod(A,listMethods=c("nipals","mean"),patternNA=rep(0.1,2),seed=1:20,nDatasets = 1)

test_that("whichNAmethod_1",{expect_true(all.equal(resultComparison,resultComparison2))})


# test on Russett
data(Russett)
X_agric =as.matrix(Russett[,c("gini","farm","rent")])
X_ind = as.matrix(Russett[,c("gnpr","labo")])
X_polit = as.matrix(Russett[ , c("demostab", "dictator")])
A = list(X_agric, X_ind, X_polit)

# Test on russett with only one variable
data(Russett)
X_agric =as.matrix(Russett[,c("gini","farm","rent")])
X_ind = as.matrix(Russett[,c("gnpr","labo")])
X_polit = as.matrix(Russett[ , c("demostab")])
A = list(agri=X_agric, ind=X_ind,polit= X_polit)
A_miss=createNA(A)
determine_patternNA(A_miss$dat)

#resultComparison2=whichNAmethod(A,rgccaType="c",listMethods=c("complete","mean","nipals","singleBlock","allInclusive"),typeNA="ponc",patternNA=rep(0.1,3),seed=1,nDatasets = 20)
#plotWhichNAmethod(listFinale=resultComparison2,ylim=c(0,0.2),output="a")

#resultComparison2=whichNAmethod(A,rgccaType="usual",typeNA="ponc",listMethods=c("complete","mean","nipals"),patternNA=rep(0.1,3),seed=NULL,nDatasets = 10)

