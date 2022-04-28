#' # print_call
#'''
data(Russett)
X_agric <- as.matrix(Russett[, c("gini", "farm", "rent")])
X_ind <- as.matrix(Russett[, c("gnpr", "labo")])
X_polit <- as.matrix(Russett[, c("demostab", "dictator")])
blocks <- list(X_agric, X_ind, X_polit)

test_that("print_call outputs the expected message", {
  res <- rgcca(blocks,
    ncomp = 1, tau = 1, scheme = "horst",
    connection = 1 - diag(3), scale = FALSE,
    scale_block = "lambda1", superblock = FALSE,
    response = NULL, NA_method = "complete"
  )
  out <- capture.output(print_call(res$call))
  expect_equal(out[1], paste0(
    "Call: method='rgcca', superblock=FALSE, scale=FALSE, ",
    "scale_block='lambda1', init='svd', bias=TRUE, tol=1e-08, ",
    "NA_method='complete', ncomp=c(1,1,1), response=NULL "
  ))
  expect_equal(out[2], "There are J = 3 blocks.")
  expect_equal(out[3], "The design matrix is:")
  expect_equal(out[4], "       block1 block2 block3")
  expect_equal(out[5], "block1      0      1      1")
  expect_equal(out[6], "block2      1      0      1")
  expect_equal(out[7], "block3      1      1      0")
  expect_equal(out[8], "")
  expect_equal(out[9], "The horst scheme is used.")
  expect_equal(length(out), 9)

  res <- rgcca(blocks,
    ncomp = 2, tau = 1, scheme = "centroid", scale = TRUE,
    scale_block = TRUE, superblock = FALSE,
    response = 3, NA_method = "nipals"
  )
  out <- capture.output(print_call(res$call))
  expect_equal(out[1], paste0(
    "Call: method='rgcca', superblock=FALSE, scale=TRUE, ",
    "scale_block=TRUE, init='svd', bias=TRUE, tol=1e-08, ",
    "NA_method='nipals', ncomp=c(2,2,2), response=3 "
  ))
  expect_equal(out[2], "There are J = 3 blocks.")
  expect_equal(out[3], "The design matrix is:")
  expect_equal(out[4], "       block1 block2 block3")
  expect_equal(out[5], "block1      0      0      1")
  expect_equal(out[6], "block2      0      0      1")
  expect_equal(out[7], "block3      1      1      0")
  expect_equal(out[8], "")
  expect_equal(out[9], "The centroid scheme is used.")
  expect_equal(length(out), 9)
})
