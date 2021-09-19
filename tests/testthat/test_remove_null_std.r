#'# remove_null_std test

#'''
 df = sapply(seq(3), function(x) runif(10))
 df = cbind(df, rep(1, 10))
 vec=rep(1,10);vec[1:4]=NA
 df3=cbind(df,vec)
 df2=remove_null_sd(list(df))$list_m
 df4=remove_null_sd(list(df3))$list_m

 test_that("remove_null_std",{expect_true(dim(df2[[1]])[2]==3)})
 test_that("remove_null_std_with_na",{expect_true(dim(df4[[1]])[2]==3)})
