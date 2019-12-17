#'# remove_null_std test

#'''
 df = sapply(seq(3), function(x) runif(10))
 df = cbind(df, rep(1, 10))
 df2=remove_null_sd(list(df))
 test_that("remove_null_std",{expect_true(dim(df2[[1]])[2]==3)})
 