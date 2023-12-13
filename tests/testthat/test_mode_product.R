test_that("mode_product works the same way as the matrix product", {
  X <- matrix(rnorm(30 * 12), nrow = 30, ncol = 12)
  y <- rnorm(12)
  res <- mode_product(X, y, m = 2)
  expect_equal(length(res), 30)
  expect_equal(res, X %*% y)

  y <- matrix(rnorm(12 * 17), nrow = 12, ncol = 17)
  res <- mode_product(X, y, m = 2)
  expect_equal(dim(res), c(30, 17))
  expect_equal(res, X %*% y)
})

test_that("mode_product works between a tensor and a vector", {
  X <- array(rnorm(40 * 5 * 7), dim = c(40, 5, 7))
  y <- rnorm(5)
  res <- mode_product(X, y, m = 2)
  expect_equal(dim(res), c(40, 1, 7))
  expect_equal(
    as.vector(res), drop(t(as.vector(X)) %*% (diag(7) %x% y %x% diag(40)))
  )

  y <- rnorm(7)
  res <- mode_product(X, y, m = 3)
  expect_equal(dim(res), c(40, 5, 1))
  expect_equal(
    as.vector(res), drop(t(as.vector(X)) %*% (y %x% diag(5) %x% diag(40)))
  )
})

test_that("mode_product works between a tensor and a matrix", {
  X <- array(rnorm(20 * 6 * 4), dim = c(20, 6, 4))
  y <- matrix(rnorm(6 * 7), nrow = 6, ncol = 7)
  res <- mode_product(X, y, m = 2)
  expect_equal(dim(res), c(20, 7, 4))
  expect_equal(
    matrix(aperm(res, c(2, 1, 3)), nrow = 7),
    t(y) %*% matrix(aperm(X, c(2, 1, 3)), nrow = 6)
  )

  y <- matrix(rnorm(4 * 13), nrow = 4, ncol = 13)
  res <- mode_product(X, y, m = 3)
  expect_equal(dim(res), c(20, 6, 13))
  expect_equal(
    matrix(aperm(res, c(3, 1, 2)), nrow = 13),
    t(y) %*% matrix(aperm(X, c(3, 1, 2)), nrow = 4)
  )
})
