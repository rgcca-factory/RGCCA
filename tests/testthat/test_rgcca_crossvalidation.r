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
        rgcca_cv_k(rgcca_res=rgcca_out,n_cores=1,validation="loo"),
         0.495311)
 }
 )
 
blocks2=blocks
blocks2[[1]][,1]=0
blocks2[[1]][1,1]=1
rgcca_out2= rgcca(blocks2, response = 1,superblock=FALSE,ncomp=1,scale=TRUE,scale_block=TRUE)
v_inds <- sample(nrow(rgcca_out2$call$raw[[1]]))
v_inds <- split(v_inds, sort(v_inds %% 5))
newA = lapply(rgcca_out2$call$raw, function(x) x[v_inds[[1]], , drop = FALSE])
# res_predict = rgcca_predict(rgcca_out2, newA = newA, model = "classification", fit = "lda", bloc_to_pred = names(rgcca_out2$call$blocks)[rgcca_out2$call$response])
# rescv= rgcca_cv_k(rgcca_res=rgcca_out2,n_cores=1,validation="loo")
 #new_scaled = FALSE si les blocs en entrée de newA ne sont pas scalés, TRUE si les blocks sont scalés
 # Finding back 0.495311
 res=rep(NA,dim(blocks[[1]])[1])
 for(i in 1:dim(blocks[[1]])[1])
 {
     # Etape 1: on extrait la ligne i des blocs non scalés
     A_moins_i=lapply(blocks,function(x){return(x[-i,])})
     A_i=lapply(blocks,function(x){return(x[i,])})
     names(A_moins_i)=names(A_i)=names(blocks)
     # on calcule la RGCCA sur le bloc A sans le i
     rgcca_out_i <- rgcca(A_moins_i, response = 1,superblock=FALSE,ncomp=1,scale=TRUE,scale_block=TRUE)
     
     # Methode 1 on predit la valeur i à l'aide du modèle rgcca, et de la ligne du bloc
     # Ici, on reprend les scale value de la rgcca calculée sur A_moins_i et on les applique à A_i
     respred_i=rgcca_predict(rgcca_out_i, A_i,new_scaled=FALSE,bloc_to_pred="agriculture") 
     
     # Methode 2: on predit la valeur i a l'aide des blocs "scalés" 
    # newA_i = lapply(rgcca_out$call$blocks, function(x) x[inds, , drop = FALSE])
    # newA_i=(rgcca_out$call$blocks)*()+() -()
   #  respred_i2=rgcca_predict(rgcca_out_i,newA_i ,new_scaled=TRUE,bloc_to_pred="agriculture") 
     
    # all.equal(respred_i,respred_i2)
    # A_i_scaled=lapply(rgcca_out$call$blocks,function(x){res=t(as.matrix(x[i,]));rownames(res)=rownames(rgcca_out$call$blocks)[1];return(res)})

    

     
       #= newA_i=scaled_A[i,]
     
     
         res[i]=respred_i$score 
    }
 mean(res)

 
 
 
# Test du predict du premier individu en leave-one-out
#--------------------------------------------------------
 rgcca_out <- rgcca(blocks, response = 1,superblock=FALSE,ncomp=1,scale=TRUE,scale_block=TRUE)
 rgcca_cv=rgcca_cv_k(rgcca_res=rgcca_out,n_cores=1,validation="loo",tol=1e-5)
 rgcca_cv$list_pred[[1]]
 # test for the first neighbor in rgcca
 rgcca_k <- set_rgcca(rgcca_out, inds = 1,tol=1e-5)
 
  newA=lapply(blocks,function(x){return(x[1,,drop=FALSE])})
 res_pred=rgcca_predict(rgcca_k,newA=newA,bloc_to_pred = "agriculture", new_scaled = FALSE)
 res_pred$prediction
 res_pred_score=res_pred$score
 test_that("rgcca_cv_k_rmse",{expect_true(
     round(res_pred_score,digits=5)== round(rgcca_cv$list_scores[1],digits=5)
 )})
 

 # Cross-validation to find out the prediction error when agri is response and with 2 comp
set.seed(1)
 rgcca_out <- rgcca(blocks, response = 1,superblock=FALSE,ncomp=2,scale=TRUE,scale_block=TRUE)
 test_that("rgcca_cv_default_2", {
          test_structure_cv(
             rgcca_cv_k(rgcca_res=rgcca_out,n_cores=1,validation="loo"),
             0.4851892)
 }
     
 )

 # k-fold
 #----------
 rgcca_out <- rgcca(blocks, response = 1,ncomp=1)

 test_that("rgcca_cv_with_kfold", {
     test_structure_cv(
         rgcca_cv_k(
             rgcca_res=rgcca_out,
         validation = "kfold",
             k = 5,
             n_cores = 1),
         scores=NA,
         val=FALSE)
 }
)

 t0= Sys.time()
 rgcca_cv_k(
     rgcca_res=rgcca_out,
     validation = "kfold",
     k = 5,
     n_cores = 1,parallelization=FALSE)
 Sys.time()-t0
 
 t0= Sys.time()
 rgcca_cv_k(
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
 #          rgcca_cv_k(
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
     cv = rgcca_cv_k(rgcca_res=rgcca_out, n_cores = 1)
     
   #  rgcca_out <-rgcca(blocksNA, response=1,ncomp=1,method="mean")
   #  cv = rgcca_cv_k(rgcca_res=rgcca_out, n_cores = 1)
     
     rgcca_out <-rgcca(blocksNA, response=1,ncomp=1,method="complete")
     cv = rgcca_cv_k(rgcca_res=rgcca_out, n_cores = 1)
     # avec la method complete -> ne fonctionne pas #TODO ?
  
 #    rgcca_out <- rgcca(blocksNA, response = 1, tol = 1E-03,method="nipals")
     cv = rgcca_cv_k(rgcca_out, n_cores = 1,validation="loo")
     test_that("rgcca_cv_withNA2", {
         round(cv$scores,digits=3)==0.443
     })
    
# Test avec scale=FALSE, scale_block=FALSE
     rgcca_out_1 <- rgcca(blocks, response = 1,superblock=FALSE,ncomp=1,scale=TRUE,scale_block=TRUE)
     rescv1=rgcca_cv_k(rgcca_res=rgcca_out_1,n_cores=1,validation="loo")
     
 
     blocks2=RGCCA:::scaling(blocks,scale=TRUE,scale_block=TRUE)
     rgcca_out_2 <- rgcca(blocks2, response = 1,superblock=FALSE,ncomp=1,scale=TRUE,scale_block=TRUE)
     rescv2=rgcca_cv_k(rgcca_res=rgcca_out_2,n_cores=1,validation="loo",scale_block=FALSE,scale=FALSE)
     
     all.equal(rescv1,rescv2)
     
     
     # rgcca_cv_k for classification
     #=====================================
     blocks_for_classif = list(
         agriculture = Russett[, 1:3],
         industry = Russett[, 4:5],
         politic = matrix(Russett[, 11],ncol=1)
     )
     blocks_for_classif[["politic"]][blocks_for_classif[["politic"]][,1]==1,]="demo"
     blocks_for_classif[["politic"]][blocks_for_classif[["politic"]][,1]==0,]="ndemo"
     
     
     
     A=blocks_for_classif
     object1 = rgcca(A, connection = C, tau = c(1,1,1),
                     ncomp = 1, superblock = FALSE, response = 3)
     rescv1=rgcca_cv_k(rgcca_res=object1,n_cores=1,validation="loo",model="classification",fit="lda")
     #   res_test  = rgcca_predict(object1, A_test,new_scaled=FALSE,fit="lda",model="classification",bloc_to_pred="politic") 
     test_that("rgcca_predict_classif",{expect_true(
         round(rescv1$score,digits=3)==0.213
     )})
     
     
     res=rep(NA,dim( blocks_for_classif[[1]])[1])
     for(i in 1:dim( blocks_for_classif[[1]])[1])
     {
         # Etape 1: on extrait la ligne i des blocs non scalés
         A_moins_i=lapply( blocks_for_classif,function(x){return(x[-i,])})
         A_i=lapply( blocks_for_classif,function(x){return(x[i,])})
         names(A_moins_i)=names(A_i)=names( blocks_for_classif)
         # on calcule la RGCCA sur le bloc A sans le i
         rgcca_out_i <- rgcca(A_moins_i, response = 3,superblock=FALSE,ncomp=1,scale=TRUE,scale_block=TRUE)
         respred_i=rgcca_predict(rgcca_out_i, A_i,new_scaled=FALSE,bloc_to_pred="politic",model="classification",fit="lda") 
         res[i]=respred_i$score 
     }
     mean(res) # 0.213