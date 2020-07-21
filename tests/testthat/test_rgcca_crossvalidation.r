set.seed(1)

data("Russett")
blocks <- list(
    agriculture = Russett[, seq(3)],
    industry = Russett[, 4:5],
    politic = Russett[, 6:11])



 test_structure_cv <- function(res, scores, nrow = 47,val=TRUE){
     expect_equal(length(res), 8)
     expect_is(res, "cv")
     expect_is(res$rgcca, "rgcca")
     pred <- res$preds
     expect_is(pred, "list")
     expect_is(pred[[1]], "matrix")
     #expect_true(all(sapply(pred, NCOL) == 2))
     expect_true(all(sapply(pred, NROW) == nrow))
      expect_identical(res$rgcca,rgcca_out)
      if(val)
          expect_identical(round(res$scores, 4), round(scores,4))
 }
 
# Crossvalidation with leave-one-out 
#----------------------------------------
 # Cross-validation to find out the prediction error when agri is response and with 1 comp
 rgcca_out <- rgcca(blocks, response = 1,superblock=FALSE,ncomp=1,scale=TRUE,scale_block=TRUE)
 test_that("rgcca_cv_default_1", {
     test_structure_cv(
         rgcca_crossvalidation(rgcca_res=rgcca_out,n_cores=1,validation="loo"),
         0.495311)
 }
 )
 
# Test du predict du premier individu en leave-one-out
#--------------------------------------------------------
 rgcca_out <- rgcca(blocks, response = 1,superblock=FALSE,ncomp=1,scale=TRUE,scale_block=TRUE)
 rgcca_cv=rgcca_crossvalidation(rgcca_res=rgcca_out,n_cores=1,validation="loo",tol=1e-5)
 rgcca_cv$list_scores
 # test for the first neighbor in rgcca
 rgcca_k <-
     RGCCA:::set_rgcca(rgcca_out,
               inds = 1,tol=1e-5)
  newA=lapply(blocks,function(x){return(x[1,,drop=FALSE])})
 res_pred=rgcca_predict(rgcca_k,newA=newA,bloc_to_pred = "agriculture", new_scaled = FALSE)$score
 
 test_that("rgcca_crossvalidation_rmse",{expect_true(
     round(res_pred,digits=5)== round(rgcca_cv$list_scores[1],digits=5)
 )})
 

 
 res=rep(NA,dim(blocks[[1]])[1])
 for(i in 1:dim(blocks[[1]])[1])
 {
     A_moins_i=lapply(blocks,function(x){return(x[-i,])})
     A_i=lapply(blocks,function(x){return(x[i,])})
     names(A_moins_i)=names(A_i)=names(blocks)
     rgcca_out_i <- rgcca(A_moins_i, response = 1,superblock=FALSE,ncomp=1,scale=TRUE,scale_block=TRUE)
     res[i]=rgcca_predict(rgcca_out_i,A_i ,new_scaled=FALSE,bloc_to_pred="agriculture") $score
 }
round(res-rgcca_cv$list_scores,digits=4)
 
 
 # Cross-validation to find out the prediction error when agri is response and with 2 comp
set.seed(1)
 rgcca_out <- rgcca(blocks, response = 1,superblock=FALSE,ncomp=2,scale=TRUE,scale_block=TRUE)
 test_that("rgcca_cv_default_2", {
          test_structure_cv(
             rgcca_crossvalidation(rgcca_res=rgcca_out,n_cores=1,validation="loo"),
             0.4851892)
 }
     
 )

 # k-fold
 #----------
 rgcca_out <- rgcca(blocks, response = 1,ncomp=1)

 test_that("rgcca_cv_with_kfold", {
     test_structure_cv(
         rgcca_crossvalidation(
             rgcca_res=rgcca_out,
         validation = "kfold",
             k = 5,
             n_cores = 1),
         scores=NA,
         val=FALSE)
 }
)

 t0= Sys.time()
 rgcca_crossvalidation(
     rgcca_res=rgcca_out,
     validation = "kfold",
     k = 5,
     n_cores = 1,parallelization=FALSE)
 Sys.time()-t0
 
 t0= Sys.time()
 rgcca_crossvalidation(
     rgcca_res=rgcca_out,
     validation = "kfold",
     n_cores=1,
     k = 5)
 Sys.time()-t0
 
# test / train #TODO
#---------------------
 # rgcca_out <- rgcca(blocks, response = 1,ncomp=1)
 # test_that("rgcca_cv_test", {
 #     test_structure_cv(
 #          rgcca_crossvalidation(
 #             rgcca_res=rgcca_out, 
 #              validation = "test", 
 #              n_cores = 1), 
 #        scores=NA,
 #        val=FALSE
 #        ) # TODO : warnings
 #     }
 # )
# 

 # Cross-validation with missing values
 #----------------------------------------
RussettWithNA <- Russett
     RussettWithNA[1:2,1:3] <- NA
     RussettWithNA[3,4:5] <- NA
     RussettWithNA[3,1] <- NA
     blocksNA <- list(
         agriculture = RussettWithNA[, seq(3)], 
         industry = RussettWithNA[, 4:5],
         politic = RussettWithNA[, 6:11])

   #  rgcca_out <- rgcca(blocksNA, response = 1,ncomp=1)
     rgcca_out <-rgcca(blocksNA, response=1,ncomp=1,method="nipals")
     cv = rgcca_crossvalidation(rgcca_res=rgcca_out, n_cores = 1)
     
     rgcca_out <-rgcca(blocksNA, response=1,ncomp=1,method="mean")
     cv = rgcca_crossvalidation(rgcca_res=rgcca_out, n_cores = 1)
     
     rgcca_out <-rgcca(blocksNA, response=1,ncomp=1,method="complete")
     cv = rgcca_crossvalidation(rgcca_res=rgcca_out, n_cores = 1)
     
     # avec la method complete -> ne fonctionne pas
  
 #    rgcca_out <- rgcca(blocksNA, response = 1, tol = 1E-03,method="nipals")
     cv = rgcca_crossvalidation(rgcca_out, n_cores = 1)
 #    test_that("rgcca_cv_withNA2", {
 #        round(cv$scores,digits=3)==0.481
 #    })
    
