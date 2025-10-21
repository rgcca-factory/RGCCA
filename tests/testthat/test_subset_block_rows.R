# '# Test subset_block_rows
#
# '''
test_that("subset_block_rows successfully extract rows of vectors", {
  x <- 1:10
  expect_equal(subset_block_rows(x, c(3, 5, 6)), x[c(3, 5, 6)])
  names(x) <- paste0("V", 1:10)
  expect_equal(subset_block_rows(x, c(3, 5, 6)), x[c(3, 5, 6)])
})

test_that("subset_block_rows successfully extract rows of matrices", {
  x <- matrix(1:21, 7, 3)
  expect_equal(subset_block_rows(x, c(3, 5, 6)), x[c(3, 5, 6), ])
  rownames(x) <- paste0("R", 1:7)
  colnames(x) <- paste0("C", 1:3)
  expect_equal(subset_block_rows(x, c(3, 5, 6)), x[c(3, 5, 6), ])
})

test_that("subset_block_rows successfully extract rows of arrays", {
  x <- array(1:72, dim = c(6, 3, 2, 2))
  expect_equal(subset_block_rows(x, c(3, 5, 6)), x[c(3, 5, 6), , , ])
  dimnames(x)[[1]] <- paste0("A", 1:6)
  dimnames(x)[[2]] <- paste0("B", 1:3)
  dimnames(x)[[3]] <- paste0("C", 1:2)
  dimnames(x)[[4]] <- paste0("D", 1:2)
  expect_equal(subset_block_rows(x, c(3, 5, 6)), x[c(3, 5, 6), , , ])
})

test_that("subset_block_rows successfully extract rows of data.frames", {
  x <- as.data.frame(matrix(1:21, 7, 3))
  expect_equal(subset_block_rows(x, c(3, 5, 6)), x[c(3, 5, 6), ])
})
