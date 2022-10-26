#' # remove_null_sd test

#'''
X1 <- vapply(seq(3), function(x) rnorm(10), FUN.VALUE = double(10))
X2 <- rep(1, 10)
X3 <- rep(1, 10)
X3[1:4] <- NA
X4 <- rep(NA, 10)
df <- cbind(X1, X2, X3)
df2 <- cbind(X1, X2, X4)

test_that("remove_null_sd removes columns with null variance", {
  res <- remove_null_sd(list(df))
  expect_equal(df[, seq(3)], res$list_m[[1]])
  removed_columns <- c(4, 5)
  names(removed_columns) <- c("X2", "X3")
  expect_equal(removed_columns, res$column_sd_null[[1]])
})

test_that("remove_null_sd removes columns of NA", {
  res <- remove_null_sd(list(df2))
  expect_equal(df[, seq(3)], res$list_m[[1]])
  removed_columns <- c(4, 5)
  names(removed_columns) <- c("X2", "X4")
  expect_equal(removed_columns, res$column_sd_null[[1]])
})
