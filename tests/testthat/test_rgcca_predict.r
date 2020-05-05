set.seed(1)
# uilding the blocks
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
 A = lapply(blocks, function(x) x[1:32,]); A[[1]][1,]
 attributes(scale(A[[1]]))$'scaled:center'
 attributes(scale(A[[1]]))$'scaled:scale'
 A_scaled=lapply(A,scale); A_scaled[[1]][1,]
 
 (A_scaled[[1]]/sqrt(3))[1,]
 
res2= scl_fun(A[[1]], attributes(scale(A[[1]]))$'scaled:center', attributes(scale(A[[1]]))$'scaled:scale')
 object1 = rgcca(A, connection = C, tau = c(0.7,0.8,0.7),
     ncomp = c(3,2,4), superblock = FALSE, response = 3)
 res2[1,]
 res  = rgcca_predict(object1, A,new_scaled=FALSE) 

test_that("rgcca_predict",{expect_true(sum(!abs(res$pred[[1]]- object1$Y[[1]])<1e-12)==0)})

 
 
 
 
 newA = lapply(blocks, function(x) x[-c(1:32),])
 newA = lapply( newA, function(x) x[, sample(1:NCOL(x))] )
 newA = sample(newA, length(newA))

 res  = rgcca_predict(object1, newA) 
 names(res)
 res$pred[[1]]
 
 mean1=apply(A[[names(blocks)[1]]],2,mean)
 sd1=apply(A[[names(blocks)[1]]],2,sd)

  mean1+sd1*(as.matrix(scale(newA[[names(blocks)[1]]]))%*%object1$a[[1]])[15,]
  res$pred[[2]][15,]
 bloc_to_pred = "industry"
