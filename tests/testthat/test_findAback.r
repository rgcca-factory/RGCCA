
# test on random data
#---------------------
set.seed(42);X1=matrix(rnorm(350),70,5);X2=matrix(rnorm(280),70,4);colnames(X1)=paste("A",1:5);colnames(X2)=paste("B",1:4);
rownames(X1)=rownames(X2)=paste("S",1:70)
A=list(X1,X2);
Ascaled=lapply(A,scale2)
B=findAback(Ascaled,sameBlockWeight=FALSE)
all.equal(A,B)

set.seed(42);X1=matrix(rnorm(350),70,5);X2=matrix(rnorm(280),70,4);colnames(X1)=paste("A",1:5);colnames(X2)=paste("B",1:4);
rownames(X1)=rownames(X2)=paste("S",1:70)
X1[1,1]=NA
X1[2,]=NA
X2[1,]=NA
X2[10:25,]=NA
Ana=list(X1,X2);
Anascaled=lapply(Ana,scale2)
B=findAback(Anascaled,sameBlockWeight=FALSE)
all.equal(Ana,B)

set.seed(42);X1=matrix(rnorm(350),70,5);X2=matrix(rnorm(280),70,4);colnames(X1)=paste("A",1:5);colnames(X2)=paste("B",1:4);
rownames(X1)=rownames(X2)=paste("S",1:70)
A=list(X1,X2);
Ascaled=lapply(A,scale2)
Ascaled2=lapply(A,function(x)
        {
            y=x/sqrt(NCOL(x));
            attr(y,"scaled:scale")=attributes(x)$'scaled:scale';
            attr(y,"scaled:center")=attributes(x)$'scaled:center' 
            return(y)
        })
B=findAback(Ascaled2)
all.equal(Ana,B)

#test_that("test_findAback",{expect_true(all.equal(A,resA)})
