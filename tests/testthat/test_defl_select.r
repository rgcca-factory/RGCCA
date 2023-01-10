### Test defl_select
set.seed(42)
X1 <- matrix(rnorm(15), 3, 5)
X2 <- matrix(rnorm(12), 3, 4)
X3 <- matrix(rnorm(21), 3, 7)
A <- list(X1, X2, X3)
yy <- cbind(c(1, 0, 0), c(0, 0, 1), c(1 / sqrt(2), 1 / sqrt(2), 0))
yy <- lapply(seq_len(NCOL(yy)), function(i) yy[, i])

test_that("defl_select does not deflate response block", {
  res <- defl_select(
    yy = yy, rr = A, nncomp = c(1, 1, 1), nn = 1, nbloc = 3, response = 3
  )
  expect_equal(res$resdefl[[3]], A[[3]])
})

test_that("defl_select does not deflate block which reached ncomp", {
  res <- defl_select(
    yy = yy, rr = A, nncomp = c(1, 2, 2), nn = 2, nbloc = 3
  )
  expect_equal(res$resdefl[[1]], A[[1]])
})

test_that("defl_select outputs coherent residuals and projections", {
  res <- defl_select(yy = yy, rr = A, nncomp = c(1, 1, 1), nn = 1, nbloc = 3)
  for (j in seq_along(A)) {
    expect_equal(A[[j]] - yy[[j]] %*% t(res$pdefl[[j]]), res$resdefl[[j]])
  }
})
