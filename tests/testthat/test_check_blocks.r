#'# check_blocks test

#'''
#  names of blocks
 rand_mat <- function(x) matrix(runif(9), 3, 3)
 A = lapply(1:3, rand_mat)
 A[[4]]=c(1,2,-1)
 try(check_blocks(A))
 geterrmessage()
 names(A) <- LETTERS[1:3]
 # importance of colnames
 try(check_blocks(A))
 geterrmessage()
 A=lapply(A,as.matrix)
 for(i in 1:length(A))
 {
     colnames(A[[i]]) <- paste0("Block",i,"_var_",letters[1:(dim(A[[i]])[2])])
 }

 # Rownames should be present   
A2= check_blocks(A)
check_blocks(A2) 

# when the blocks have not the same dimensions
X1=matrix(rnorm(15),5,3);rownames(X1)=paste("S",1:5);colnames(X1)=letters[1:3]
X2=matrix(rnorm(16),4,4);rownames(X2)=paste("S",c(1,2,5,4));colnames(X2)=letters[5:8]
X3=matrix(rnorm(3),3,1);rownames(X3)=paste("S",c(5,1,3));colnames(X3)="z"
A=list(X1,X2,X3)
A2=check_blocks(A)
A2=try(check_blocks(A,add_NAlines = TRUE))
test_that("check_blocks",{expect_true(A2[[3]]["S 5",1]==A[[3]]["S 5",1])})


 A[[1]][2, 3] <- "character"
 try(check_blocks(A))
# A[[1]][2, 3] <- runif(1)
# init : boolean (FALSE by default) for the first block checking


# test with blocks with 1 column

