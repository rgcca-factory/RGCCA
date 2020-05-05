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
 A = lapply(blocks, function(x) x[1:32,]);
 
object1 = rgcca(A, connection = C, tau = c(0.7,0.8,0.7),
     ncomp = c(3,2,4), superblock = FALSE, response = 3)
 res  = rgcca_predict(object1, A,new_scaled=FALSE) 
 
# Testing res$pred
test_that("rgcca_predict",{expect_true(sum(!abs(res$pred[[1]]- object1$Y[[1]])<1e-12)==0)})

 
 