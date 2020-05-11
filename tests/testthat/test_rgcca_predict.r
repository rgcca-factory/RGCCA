set.seed(1)
# Building the blocks
 data("Russett")
 blocks = list(
 agriculture = Russett[, 1:3],
 industry = Russett[, 4:5],
 politic = Russett[, 6:11]
 )
 C = connection = matrix(c(0, 0, 1,
 0, 0, 1,
 1, 1, 0),
 3, 3)
 A = lapply(blocks, function(x) x[1:32,]);
 newA=lapply(A,scale)
 (newA[[1]]/sqrt(3))[1,]

object1 = rgcca(A, connection = C, tau = c(0.7,0.8,0.7),
     ncomp = c(3,2,4), superblock = FALSE, response = 3)
apply(object1$call$block[[1]],2,sd)
 res  = rgcca_predict(object1, A,new_scaled=FALSE) 
test_that("rgcca_predict",{expect_true(
    sum(!abs(res$pred[[1]]- object1$Y[[1]])<1e-12)==0
    )})


A_restr=lapply(blocks,function(x) x[1:16,])
res_restr  = rgcca_predict(object1, A_restr,new_scaled=FALSE) 
test_that("rgcca_predict_restr",{expect_true(
    sum(!abs(res_restr$pred[[1]]["Argentina",]- object1$Y[[1]]["Argentina",])<1e-12)==0
    )})

# Testing on a different dataset
A_test=lapply(blocks,function(x) x[c(1,33),])
res_test  = rgcca_predict(object1, A_test,new_scaled=FALSE) 
test_that("rgcca_predict_two",{expect_true(
        sum(!abs(res_test$pred[[1]]["Argentina",]- object1$Y[[1]]["Argentina",])<1e-12)==0
        )})
test_that("rgcca_predict_two2",{expect_true(round(res_test$pred[[1]]["Panama",1],digits=4)==0.07)})

#With one line only
A_test=lapply(blocks,function(x) x[c(33),])
res_test_one  = rgcca_predict(object1, A_test,new_scaled=FALSE) 
test_that("rgcca_predict_one_ind",{expect_true(round(res_test_one$pred[[1]]["Panama",1],digits=4)==0.07)})

# With a block to predict
A = lapply(blocks, function(x) x[1:32,]);
newA=lapply(A,scale)
(newA[[1]]/sqrt(3))[1,]
object1 = rgcca(A, connection = C, tau = c(0.7,0.8,0.7),
                ncomp = c(3,2,4), superblock = FALSE, response = 3)
res  = rgcca_predict(object1, A,new_scaled=FALSE,bloc_to_pred="politic") 

# With a univariate block to predict
A = lapply(blocks, function(x) x[1:32,]);
newA=lapply(A,scale)
(newA[[1]]/sqrt(3))[1,]
object1 = rgcca(A, connection = C, tau = c(0.7,0.8,0.7),
                ncomp = c(1,1,1), superblock = FALSE, response = 3)
res  = rgcca_predict(object1, A,new_scaled=FALSE,bloc_to_pred="politic") 

# Checking the "ordre de grandeur"
#res$pred_A
#head(res$pred_A)
head(A[[3]])
res$res
test_that("rgcca_predict_with_block",{expect_true(sum(!abs(res$pred[[1]]- object1$Y[[1]])<1e-12)==0)})


# cor
#A = lapply(blocks, function(x) x[1:32,]);
#A_test=lapply(blocks,function(x) x[c(39:47),])
#object1 = rgcca(A, connection = C, tau = c(0.7,0.8,0.7),
#                ncomp = c(3,2,4), superblock = FALSE, response = 3)

res_test  = rgcca_predict(object1, A_test,new_scaled=FALSE,fit="cor",bloc_to_pred="politic") 
#test_that("rgcca_predict_two",{expect_true(
#    sum(!abs(res_test$pred[[1]]["Sweden",]- object1$Y[[1]]["Sweden",])<1e-12)==0
#)})

# TODO
#A_test=lapply(blocks,function(x) x[c(39:47),])
#res_test  = rgcca_predict(object1, A_test,new_scaled=FALSE,fit="lda",model="classification",bloc_to_pred="politic") 
#test_that("rgcca_predict_two",{expect_true(
#    sum(!abs(res_test$pred[[1]]["Argentina",]- object1$Y[[1]]["Argentina",])<1e-12)==0
#)})




 # With missing values
RussettWithNA <- Russett
RussettWithNA[1:2,1:3] <- NA
RussettWithNA[3,4:5] <- NA
RussettWithNA[3,1] <- NA
blocksNA <- list(
    agriculture = RussettWithNA[, seq(3)], 
    industry = RussettWithNA[, 4:5],
    politic = RussettWithNA[, 6:11])
A_test=lapply(blocksNA,function(x) x[c(39:47),])
object1 = rgcca(A, connection = C, tau = c(0.7,0.8,0.7),
                ncomp = c(1,1,1), superblock = FALSE, response = 3)

res_test  = rgcca_predict(object1, A_test,new_scaled=FALSE,bloc_to_pred="politic") 

