data(Russett)
 blocks = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
     politic = Russett[, 6:11] )
 rgcca_out = rgcca(blocks,ncomp=2)
 t=get_cor_all(rgcca_out)

 test_that("get_cor_all",
           {
               expect_true(length(t)==2)
           })
