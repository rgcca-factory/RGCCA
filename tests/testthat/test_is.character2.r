#'# is.character2 test

#'''
 x = matrix(c(runif(10), LETTERS[seq(10)]), 10, 2)
 is.character2(x)
 # TRUE
 is.character2(LETTERS[seq(10)])
 # TRUE
test_that("is.character2",{expect_true(is.character2(LETTERS[seq(10)]))})
test_that("is.character2",{expect_true( is.character2(x))})
