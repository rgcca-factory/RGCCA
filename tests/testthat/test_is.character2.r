#' # is.character2 test

#'''
x <- matrix(c(rnorm(10), LETTERS[seq(10)]), 10, 2)
test_that("is.character2 returns TRUE if x contains at least a character", {
  expect_true(is.character2(LETTERS[seq(10)]))
  expect_true(is.character2(x))
})
test_that("is.character2 returns FALSE if x does not contain any character", {
  expect_false(is.character2(1:10))
  expect_false(is.character2(c(rnorm(5), "NA", NA)))
})
