X <- crossprod(matrix(rnorm(10 * 5), 10, 5))
X_inv <- ginv(X)

test_that("sqrt_matrix retrieves the square root of a symmetric matrix", {
  sqrt_X <- sqrt_matrix(X, inv = FALSE)
  expect_equal(X, sqrt_X %*% sqrt_X)
})

test_that("sqrt_matrix retrieves the inverse of a square root of a
          symmetric matrix", {
  sqrt_inv_X <- sqrt_matrix(X, inv = TRUE)
  expect_equal(X_inv, sqrt_inv_X %*% sqrt_inv_X)
})
