#'# khatri_rao test

#'''
# Expect error when x or y is not a matrix
test_that("test_error_when_not_matrix", {
  expect_error(khatri_rao(x = matrix(1, 5, 2), y = 1:5), 
               "x and y must be matrices", 
               fixed=TRUE)
})

# Expect error when x and y do not have the same number of columns
test_that("test_error_when_different_number_of_columns", {
  expect_error(khatri_rao(x = matrix(1, 5, 2), y = matrix(1, 5, 3)), 
               "x and y must have the same number of columns", 
               fixed=TRUE)
})

# Expect to match computed by hand Khatri-Rao product
x = matrix(c(
  1,  0,  2,
  0, -1,  1
), nrow = 2, ncol = 3, byrow = T)
y = matrix(c(
   0,  3,  1,
   1,  2, -5,
  -1,  0,  6 
), nrow = 3, ncol = 3, byrow = T)
res = matrix(c(
   0,  0,  2,
   1,  0, -10,
  -1,  0,  12,
   0, -3,  1,
   0, -2, -5,
   0,  0,  6
), nrow = 6, ncol = 3, byrow = T)
test_that("test_right_result", {
  expect_identical(khatri_rao(x = x, y = y), res)
})
