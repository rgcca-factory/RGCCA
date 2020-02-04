#'# to_numeric test

#'''
X1=matrix(1:20,4,5);rownames(X1)=paste("S",1:4);colnames(X1)=paste("A",1:5)
X2=matrix(1:4,4,1);rownames(X2)=paste("S",1:4);colnames(X2)="C"
A=list(X1,X2)

# scale=TRUE, sameBlockWeight=FALSE
A1=scaling(A,bias=FALSE,scale=TRUE,sameBlockWeight=FALSE)
sd(A1[[1]][,1])==1
sd(A1[[2]][,1])==1
test_that("test_scaling_sameBWs",{expect_true(sd(A1[[1]][,1])==1)})

# scale=TRUE, sameBlockWeight=TRUE
A2a1=A1[[1]]/sqrt(5)
A2b=scaling(A,bias=FALSE,scale=TRUE,sameBlockWeight=TRUE)
test_that("test_scaling_sameBWs",{expect_true(sum(A2b[[1]]==A2a1)==20)})

# scale=FALSE, sameBlockWeight=TRUE
A3=scaling(A,bias=FALSE,scale=FALSE,sameBlockWeight=TRUE)
lapply(A3,function(x){apply(x,2,mean)})
covarMat=cov2(A[[1]],bias=FALSE); 
var(A[[1]][,1])+var(A[[1]][,2])+var(A[[1]][,3])+var(A[[1]][,4])+var(A[[1]][,5])==sum(diag(covarMat))
A3b=scale(A[[1]],scale=FALSE)/sqrt(sum(diag(covarMat)))
test_that("test_scaling_sameBW",{expect_true(sum(A3b[[1]]-A3[[1]]>1e-16)==0)})
attr(A3b[[1]],"scaled:center")
attr(A3b[[1]],"scaled:scale")

# scale=FALSE, sameBlockWeight=FALSE
A4=scaling(A,bias=FALSE,scale=FALSE,sameBlockWeight=FALSE)


A5=lapply(A,scale2)
res=scaling(A5,sameBlockWeight=TRUE,scale=FALSE,bias=FALSE)
test_that("test_scaling_end",{expect_true(sum(abs(res[[1]]-A2b[[1]])>1e-15)==0)})
 # ok
# with missing values ?
#test_that("to_numeric",{expect_true()})


