# Expect to match computed by hand Khatri-Rao product
x <- matrix(c(
  1,  0,  2,
  0, -1,  1
), nrow = 2, ncol = 3, byrow = TRUE)

y <- matrix(c(
  0,  3,  1,
  1,  2, -5,
  -1,  0,  6
), nrow = 3, ncol = 3, byrow = TRUE)

res <- matrix(c(
  0,  0,  2,
  1,  0, -10,
  -1,  0,  12,
  0, -3,  1,
  0, -2, -5,
  0,  0,  6
), nrow = 6, ncol = 3, byrow = T)

test_that("test_khatri_rao", {
  expect_identical(khatri_rao(x = x, y = y), res)
})
