# '#createNA3 test
# 
# '''
set.seed(42);X1=matrix(rnorm(500),100,5);X2=matrix(rnorm(400),100,4)
rownames(X1)=rownames(X2)=paste("S",1:dim(X1)[1],sep="")
A=list(X1,X2)
# createNA
jdd=createNA(A,typeNA="block",pNA=c(0.1,0.3))
test_that("createNA_1",{expect_true(sum(is.na(jdd$dat[[1]][,1]))==10)})
test_that("createNA_2",{expect_true(sum(is.na(jdd$dat[[2]][,1]))==30)})



jdd$dat[[1]][jdd$subjectKept,]
jdd$dat[[2]][jdd$subjectKept,]


res=createNA(A,typeNA="block",pNA=c(0.5,0.5))

jdd=createNA(A,typeNA="ponc",pNA=c(0.1,0.3))
jdd$dat[[1]][jdd$subjectKept,]
jdd$dat[[2]][jdd$subjectKept,]

res=createNA(A,typeNA="block",pNA=c(0.5,0.5),seed=30)
res2=createNA(A,typeNA="block",pNA=c(0.5,0.5),seed=30)
test_that("createNA_3",{expect_true(all.equal(res,res2))})

# with only one variable in one block
set.seed(42);X1=matrix(rnorm(500),100,5);X2=matrix(rnorm(400),100,4)
rownames(X1)=rownames(X2)=paste("S",1:dim(X1)[1],sep="")
A=list(X1,X2[,1])
# createNA
jdd=createNA(A,typeNA="block",pNA=c(0.1,0.3))

