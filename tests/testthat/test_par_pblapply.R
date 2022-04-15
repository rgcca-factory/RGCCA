test_that("par_pblapply performs lapply", {
  X <- list(rnorm(7), 1:5, c(rep(0, 3), rep(1, 9)))
  FUN <- function(x) 2 * x + x^2 - 12
  expect_identical(lapply(X, FUN), par_pblapply(X, FUN, n_cores = 1))
})
