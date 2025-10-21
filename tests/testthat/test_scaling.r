data("Russett")
blocks <- list(
  agriculture = Russett[, seq(3)],
  industry = Russett[, 4:5],
  politic = Russett[, 6:11]
)
blocks <- lapply(blocks, data.matrix)

n <- NROW(blocks[[1]])
blocks[[4]] <- array(rnorm(n * 10 * 7), dim = c(n, 10, 7))

bias <- FALSE
sqrt_N <- sqrt(n + bias - 1)
tolerance <- 1e-12

test_that("scaling_default_1", {
  blocks3 <- lapply(blocks, function(x) scale(matrix(x, nrow = nrow(x))))
  blocks3 <- lapply(blocks3, function(x) {
    return(x / sqrt(ncol(x)))
  })
  blocks3[[4]] <- array(blocks3[[4]], dim = dim(blocks[[4]]))
  blocks2 <- scaling(blocks, scale = TRUE, scale_block = TRUE, bias = bias)
  expect_true(sum(abs(blocks3[[2]] - blocks2[[2]])) < tolerance)
})

test_that("scale_block = 'inertia' leads to unit Frobenius norm", {
  b <- scaling(blocks, scale = TRUE, scale_block = TRUE, bias = bias)
  for (j in seq_along(b)) {
    expect_equal(
      norm(matrix(b[[j]], nrow = nrow(b[[j]])) / sqrt_N, type = "F"),
      1, tolerance = tolerance
    )
  }
  b <- scaling(blocks, scale = FALSE, scale_block = TRUE, bias = bias)
  for (j in seq_along(b)) {
    expect_equal(
      norm(matrix(b[[j]], nrow = nrow(b[[j]])) / sqrt_N, type = "F"),
      1, tolerance = tolerance
    )
  }
})

test_that("scale_block = 'lambda1' leads to top eigenvalue of covariance
          matrix being equal to one", {
  b <- scaling(blocks, scale = TRUE, scale_block = "lambda1", bias = bias)
  for (j in seq_along(b)) {
    expect_equal(eigen(crossprod(
      matrix(b[[j]], nrow = nrow(b[[j]])) / sqrt_N
    ))$values[1], 1, tolerance = tolerance)
  }
  b <- scaling(blocks, scale = FALSE, scale_block = "lambda1", bias = bias)
  for (j in seq_along(b)) {
    expect_equal(eigen(crossprod(
      matrix(b[[j]], nrow = nrow(b[[j]])) / sqrt_N
    ))$values[1], 1, tolerance = tolerance)
  }
})

test_that("another value of scale_block does not lead to further scaling", {
  b <- scaling(blocks, scale = TRUE, scale_block = "none", bias = bias)
  b_ref <- lapply(blocks, function(x) scale(matrix(x, nrow = nrow(x))))
  b_ref[[4]] <- array(b_ref[[4]], dim = dim(b[[4]]))
  for (j in seq_along(b)) {
    expect_lte(sum(abs(b[[j]] - b_ref[[j]])), tolerance)
  }
  b <- scaling(blocks, scale = FALSE, scale_block = "none", bias = bias)
  b_ref <- lapply(blocks, function(x) scale(
    matrix(x, nrow = nrow(x)), scale = FALSE)
  )
  b_ref[[4]] <- array(b_ref[[4]], dim = dim(b[[4]]))
  for (j in seq_along(b)) {
    expect_lte(sum(abs(b[[j]] - b_ref[[j]])), tolerance)
  }
})
