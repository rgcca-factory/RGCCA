test_that("scale2 only centers the data when scale = FALSE", {
  M <- matrix(rnorm(10 * 12), nrow = 10)
  expect_equal(scale(M, scale = FALSE), scale2(M, scale = FALSE))
})

test_that("scale2 centers and scales the data like stats::scale when
          scale = TRUE and bias = FALSE", {
  M <- matrix(rnorm(10 * 12), nrow = 10)
  expect_equal(scale(M, scale = TRUE), scale2(M, scale = TRUE, bias = FALSE))
})

test_that("scale2 centers and scales the data differently from stats::scale when
          scale = TRUE and bias = TRUE", {
  M <- matrix(rnorm(10 * 12), nrow = 10)
  M_scaled <- scale2(M, scale = TRUE, bias = TRUE)
  M_scaled <- M_scaled * sqrt(9 / 10)
  attr(M_scaled, "scaled:scale") <-
    attr(M_scaled, "scaled:scale") / sqrt(9 / 10)
  expect_equal(scale(M, scale = TRUE), M_scaled)
})
