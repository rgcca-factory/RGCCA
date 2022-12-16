### Test cov2

test_that("cov2 behaves like cov when bias = FALSE", {
  x <- matrix(rnorm(17 * 5), nrow = 17, ncol = 5)
  y <- matrix(rnorm(17 * 4), nrow = 17, ncol = 4)
  expect_equal(cov(x), cov2(x, bias = FALSE))
  expect_equal(cov(y), cov2(y, bias = FALSE))
  expect_equal(cov(x, y), cov2(x, y, bias = FALSE))
})

test_that("cov2 changes the normalization factor when bias = TRUE", {
  x <- matrix(rnorm(17 * 5), nrow = 17, ncol = 5)
  y <- matrix(rnorm(17 * 4), nrow = 17, ncol = 4)
  expect_equal(cov(x), cov2(x, bias = TRUE) * 17 / 16)
  expect_equal(cov(y), cov2(y, bias = TRUE) * 17 / 16)
  expect_equal(cov(x, y), cov2(x, y, bias = TRUE) * 17 / 16)
})
