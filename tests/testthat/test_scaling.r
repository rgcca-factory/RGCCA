data("Russett")
blocks <- list(
    agriculture = Russett[, seq(3)],
    industry = Russett[, 4:5],
    politic = Russett[, 6:11])

bias = FALSE
sqrt_N = sqrt(NROW(blocks[[1]]) + bias - 1)

blocks3 = lapply(blocks, scale)
blocks3 = lapply(blocks3, function(x) return(x/sqrt(ncol(x))))
blocks2 = scaling(blocks, scale = T, scale_block = T, bias = bias)
test_that("scaling_default_1", {
    expect_true(sum(abs(blocks3[[2]] - blocks2[[2]])) < 1e-14)})

test_that("scale_block = 'inertia' leads to unit Frobenius norm", {
  b = scaling(blocks, scale = T, scale_block = T, bias = bias)
  for (j in seq_along(b)) {
    expect_equal(norm(b[[j]] / sqrt_N, type = "F"), 1, tolerance = 1e-14)
  }
  b = scaling(blocks, scale = F, scale_block = T, bias = bias)
  for (j in seq_along(b)) {
    expect_equal(norm(b[[j]] / sqrt_N, type = "F"), 1, tolerance = 1e-14)
  }
})

test_that("scale_block = 'lambda1' leads to top eigenvalue of covariance
          matrix being equal to one", {
            b = scaling(blocks, scale = T, scale_block = "lambda1", bias = bias)
            for (j in seq_along(b)) {
              expect_equal(eigen(crossprod(b[[j]] / sqrt_N))$values[1],
                           1, tolerance = 1e-14)
            }
            b = scaling(blocks, scale = F, scale_block = "lambda1", bias = bias)
            for (j in seq_along(b)) {
              expect_equal(eigen(crossprod(b[[j]] / sqrt_N))$values[1],
                           1, tolerance = 1e-14)
            }
})

test_that("another value of scale_block does not lead to further scaling", {
  b = scaling(blocks, scale = T, scale_block = "none", bias = bias)
  b_ref = lapply(blocks, scale)
  for (j in seq_along(b)) {
    expect_true(sum(abs(b[[j]] - b_ref[[j]])) < 1e-14)
  }
  b = scaling(blocks, scale = F, scale_block = "none", bias = bias)
  b_ref = lapply(blocks, scale, center = T, scale = F)
  for (j in seq_along(b)) {
    expect_true(sum(abs(b[[j]] - b_ref[[j]])) < 1e-14)
  }
})
