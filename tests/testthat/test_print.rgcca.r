#' # print.rgcca
#'''
data(Russett)
X_agric <- as.matrix(Russett[, c("gini", "farm", "rent")])
X_ind <- as.matrix(Russett[, c("gnpr", "labo")])
X_polit <- as.matrix(Russett[, c("demostab", "dictator")])
blocks <- list(X_agric, X_ind, X_polit)

test_that("print.rgcca outputs the expected message", {
  res <- rgcca(blocks, ncomp = 1, tau = c(1, 0.7, 0.9))
  out <- capture.output(print(res))
  expect_equal(out[seq(9)], capture.output(print_call(res$call)))
  expect_equal(out[10], "Sum_{j,k} c_jk g(cov(X_ja_j, X_ka_k) = 1.1932 ")
  expect_equal(out[11], "")
  expect_equal(out[12], "The regularization parameter used for block1 is: 1")
  expect_equal(out[13], "The regularization parameter used for block2 is: 0.7")
  expect_equal(out[14], "The regularization parameter used for block3 is: 0.9")
  expect_equal(length(out), 14)

  tau <- matrix(c(1, 1, 1, 0, 0, 0), nrow = 2, byrow = TRUE)
  res <- rgcca(blocks, ncomp = 2, tau = tau)
  out <- capture.output(print(res))
  expect_equal(out[seq(9)], capture.output(print_call(res$call)))
  expect_equal(out[10], "Sum_{j,k} c_jk g(cov(X_ja_j, X_ka_k) = 1.7093 ")
  expect_equal(out[11], "")
  expect_equal(out[12], "The regularization parameters used are: ")
  expect_equal(out[13], "     block1 block2 block3")
  expect_equal(out[14], "[1,]      1      1      1")
  expect_equal(out[15], "[2,]      0      0      0")
  expect_equal(length(out), 15)

  res <- rgcca(blocks, ncomp = 1, sparsity = c(0.7, 1, 0.9), method = "sgcca")
  out <- capture.output(print(res))
  expect_equal(out[seq(9)], capture.output(print_call(res$call)))
  expect_equal(out[10], "Sum_{j,k} c_jk g(cov(X_ja_j, X_ka_k) = 0.9348 ")
  expect_equal(out[11], "")
  expect_equal(out[12], paste0(
    "The sparsity parameter used for block1 is: ",
    "0.7 (with 2 variables selected)"
  ))
  expect_equal(out[13], paste0(
    "The sparsity parameter used for block2 is: ",
    "1 (with 2 variables selected)"
  ))
  expect_equal(out[14], paste0(
    "The sparsity parameter used for block3 is: ",
    "0.9 (with 2 variables selected)"
  ))
  expect_equal(length(out), 14)

  sparsity <- matrix(c(0.8, 1, 0.9, 1, 0.9, 0.8), nrow = 2, byrow = TRUE)
  res <- rgcca(blocks, ncomp = 2, sparsity = sparsity, method = "sgcca")
  out <- capture.output(print(res))
  expect_equal(out[seq(9)], capture.output(print_call(res$call)))
  expect_equal(out[10], "Sum_{j,k} c_jk g(cov(X_ja_j, X_ka_k) = 0.9958 ")
  expect_equal(out[11], "")
  expect_equal(out[12], "The sparsity parameters used are: ")
  expect_equal(out[13], "     block1 block2 block3")
  expect_equal(out[14], "[1,]    0.8    1.0    0.9")
  expect_equal(out[15], "[2,]    1.0    0.9    0.8")
  expect_equal(out[16], "The number of selected variables are: ")
  expect_equal(out[17], "     block1 block2 block3")
  expect_equal(out[18], "[1,]      2      2      2")
  expect_equal(out[19], "[2,]      3      2      2")
  expect_equal(length(out), 19)
})
