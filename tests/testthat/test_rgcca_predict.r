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
                           1, 1, 0), 3, 3)

 A = lapply(blocks, function(x) x[1:32,]);

#-------------------------------------------------------------------------
# Checking the Y with the prediction with the response block in last position
#-------------------------------------------------------------------------
# With a block to predict
A = lapply(blocks, function(x) x[1:32,]);
object1 = rgcca(A, connection = C, tau = c(0.7,0.8,0.7),
                ncomp = c(3,2,4), superblock = FALSE, response = 3)
res  = rgcca_predict(rgcca_res = object1, X = A, X_scaled = FALSE, block_to_predict = "politic")
reslm_1_res=apply(object1$call$raw$politic,2,function(x){lm(x~object1$Y[[1]][,1:3]+object1$Y[[2]][,1:2])$residuals})
rmse_test=apply(reslm_1_res,2,function(x){return(sqrt(mean(x^2,na.rm=T)))})
test_that("rgcca_predict_rmse2",{expect_true(
    round(res$score,digits=5)==round(mean(rmse_test),digits=5)
)})

# Checking the RMSE criterion and the linear models
#-----------------------------------------------------
# for same dataset
A3 = lapply(blocks, function(x) x[1:32,]);
newA3=A3
object3 = rgcca(A3, connection = C, tau = c(0.7,0.8,0.7),
                ncomp = c(1,1,1), superblock = FALSE, response = 1)

res  = rgcca_predict(object3, newA3, X_scaled = FALSE, block_to_predict = "agriculture")
reslm_1_res=apply(object3$call$raw$agriculture,2,function(x){lm(x~object3$Y[[2]][,1]+object3$Y[[3]][,1])$residuals})

test_that("rgcca_predict_rmse",{expect_true(
    sum(round(res$res,digits=5)-round(reslm_1_res,digits=5))==0
)})

rmse_test=apply(reslm_1_res,2,function(x){return(sqrt(mean(x^2,na.rm=T)))})
test_that("rgcca_predict_rmse2",{expect_true(
    round(res$score,digits=5)==round(mean(rmse_test),digits=5)
)})

# for one line
A4 = lapply(blocks, function(x) x[1:32,]);
newA4=lapply(blocks,function(x){return(x[1,,drop=FALSE])})
object4 = rgcca(A4, connection = C, tau = c(0.7,0.8,0.7),
                ncomp = c(1,1,1), superblock = FALSE, response = 1)
res  = rgcca_predict(object4, newA4, X_scaled = FALSE, block_to_predict = "agriculture")
test_that("rgcca_predict_rmse3",{expect_true(
    sum(round(res$res,digits=5)==round(reslm_1_res[1,],digits=5))==3
)})



# cor
#A = lapply(blocks, function(x) x[1:32,]);
#A_test=lapply(blocks,function(x) x[c(39:47),])
#                ncomp = c(3,2,4), superblock = FALSE, response = 3)

#res_test  = rgcca_predict(object1, A_test,X_scaled=FALSE,prediction_model="cor",block_to_predict="politic")
#test_that("rgcca_predict_two",{expect_true(
#    sum(!abs(res_test$pred[[1]]["Sweden",]- object1$Y[[1]]["Sweden",])<1e-12)==0
#)})






#-----------------------
# With missing values
#-----------------------
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
res_test  = rgcca_predict(object1, A_test, X_scaled = FALSE, block_to_predict = "politic")

# classfication and LDA
#---------------------------------
# Blocks for classification
 blocks_for_classif = list(
     agriculture = Russett[, 1:3],
     industry = Russett[, 4:5],
     politic = matrix(Russett[, 11],ncol=1)
 )
 blocks_for_classif[["politic"]][blocks_for_classif[["politic"]][,1]==1,]="demo"
 blocks_for_classif[["politic"]][blocks_for_classif[["politic"]][,1]==0,]="ndemo"

 A = lapply(blocks_for_classif, function(x) x[1:32,]);
 A_test=lapply(blocks_for_classif,function(x) x[c(39:47),])
 object1 = rgcca(A, connection = C, tau = c(1,1,1),
                 ncomp = c(3,2,1), superblock = FALSE, response = 3)
 res_test  = rgcca_predict(object1, X = A, X_scaled = FALSE, prediction_model = "lda", task = "classification", block_to_predict = "politic")

#   res_test  = rgcca_predict(object1, A_test,X_scaled=FALSE,prediction_model="lda",task="classification",block_to_predict="politic")
 test_that("rgcca_predict_classif",{expect_true(
     res_test$score==1-0.875
     )})
#


# A = lapply(blocks_for_classif, function(x) x[1:40,]);
 #   object1 = rgcca(A, connection = C, tau = c(1,1,1),
 #                 ncomp = c(3,2,1), superblock = FALSE, response = 3)
 # res_test  = rgcca_predict(object1, X=A,X_scaled=FALSE,prediction_model="logistic",task="classification",block_to_predict="politic")

 #   res_test  = rgcca_predict(object1, A_test,X_scaled=FALSE,prediction_model="lda",task="classification",block_to_predict="politic")
 #test_that("rgcca_predict_classif",{expect_true(
 #    res_test$score==0.875
 #)})
