#' # imputeEM test

#'''
df <- matrix(runif(20), 10, 2)
r <- check_quantitative(df, "data")
test_that("test_check_quantitative", {
  expect_true(is.null(r))
})
