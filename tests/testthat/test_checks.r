# Test check_integer
test_that("an error is raised if x is not numeric", {
  expect_error(check_integer("x", "toto"),
               "x should be numeric.",
               fixed = TRUE)
})
test_that("an error is raised if x contains NA", {
  expect_error(check_integer("x", c(42, NA)),
               "x should not be NA.",
               fixed = TRUE)
})
test_that("an error is raised if type is scalar and x not of length 1", {
  expect_error(check_integer("x", c(42, 7), type = "scalar"),
               "x should be of length 1.",
               fixed = TRUE)
})
test_that("an error is raised if any element of x is a float but float is false", {
  expect_error(check_integer("x", c(1, 1.7, 2), type = "vector", float = FALSE),
               "x should be an integer.",
               fixed = TRUE)
})
test_that("an error is raised if any element of x is below min", {
  expect_error(check_integer("x", c(0, 1, 2), type = "vector", min = 1),
               "x should be higher than or equal to 1.",
               fixed = TRUE)
})
test_that("check_integer passes and return x when x is valid", {
  expect_equal(check_integer(1), 1)
  expect_equal(check_integer(c(1, 2, 3), type = "vector"), c(1, 2, 3))
  expect_equal(check_integer(1.7, float = TRUE), 1.7)
  x = matrix(1:4, 2, 2)
  rownames(x) = paste0("R", 1:2)
  colnames(x) = paste0("C", 1:2)
  expect_equal(check_integer(x, type = "matrix"), x)
  x = as.data.frame(x)
  expect_equal(check_integer(x, type = "data.frame"), x)
})
