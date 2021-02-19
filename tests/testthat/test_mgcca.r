# '# Test mgcca
# 
# '''
set.seed(0)

# Generating fake data
generate_block <- function(dim) {
  array(rnorm(prod(dim)), dim = dim)
}
X1 = generate_block(c(40, 20, 30))
X2 = generate_block(c(40, 35))
X3 = generate_block(c(40, 18, 25, 7))
A = list(X1, X2, X3);
C = matrix(1, nrow = 3, ncol = 3) - diag(3)

# Run directly mgcca or through calling rgcca. Both should run without errors.
fit.mgcca = mgcca(A, C, ncomp=rep(3, 3, 3), init = "svd", tau = c(1, 1, 1),
                  scheme = "factorial", scale = TRUE, verbose=FALSE,
                  ranks = c(3, 1, 2))
fit.rgcca = rgcca(blocks = A, connection = C, ncomp=rep(3, 3, 3), init = "svd", 
                  tau = c(1, 1, 1), scheme = "factorial", scale = TRUE, 
                  verbose=FALSE, type = "mgcca", ranks = c(3, 1, 2))

# Check number and dimensions of factors
test_that("test_length_factors", {
  expect_equal(length(fit.mgcca$factors), 3)
})
test_that("test_number_of_modes_match_number_of_factors", {
  expect_equal(length(fit.mgcca$factors[[1]]), 2)
  expect_true(is.null(fit.mgcca$factors[[2]]))
  expect_equal(length(fit.mgcca$factors[[3]]), 3)
})
test_that("test_dim_of_modes_match_dim_of_factors", {
  expect_equal(dim(fit.mgcca$factors[[1]][[1]]), c(20, 3 * 3))
  expect_equal(dim(fit.mgcca$factors[[1]][[2]]), c(30, 3 * 3))
  expect_equal(dim(fit.mgcca$factors[[3]][[1]]), c(18, 3 * 2))
  expect_equal(dim(fit.mgcca$factors[[3]][[2]]), c(25, 3 * 2))
  expect_equal(dim(fit.mgcca$factors[[3]][[3]]), c(7, 3 * 2))
})

# With SVD initialization, MGCCA is a deterministic algorithm so results should
# be the same
test_that("test_same_factors", {
  expect_true(all(unlist(sapply(1:2, function(x) 
    fit.mgcca$factors[[1]][[x]] == fit.rgcca$factors[[1]][[x]]
  ))))
  expect_true(all(unlist(sapply(1:3, function(x) 
    fit.mgcca$factors[[3]][[x]] == fit.rgcca$factors[[3]][[x]]
  ))))
})
test_that("test_same_blocks", {
  for (j in 1:length(A)) {
    expect_true(all(fit.mgcca$blocks[[j]] == fit.rgcca$blocks[[j]]))
  }
})
test_that("test_same_y", {
  for (j in 1:length(A)) {
    expect_true(all(fit.mgcca$Y[[j]] == fit.rgcca$Y[[j]]))
  }
})
test_that("test_same_a", {
  for (j in 1:length(A)) {
    expect_true(all(fit.mgcca$a[[j]] == fit.rgcca$a[[j]]))
  }
})
test_that("test_same_astar", {
  for (j in 1:length(A)) {
    expect_true(all(fit.mgcca$astar[[j]] == fit.rgcca$astar[[j]]))
  }
})
test_that("test_same_AVE", {
  for (j in 1:length(A)) {
    expect_true(all(unlist(fit.mgcca$AVE[[j]]) == unlist(fit.rgcca$AVE[[j]])))
  }
  expect_equal(fit.mgcca$AVE_outer_model, fit.rgcca$AVE_outer_model)
  expect_equal(fit.mgcca$AVE_inner_model, fit.rgcca$AVE_inner_model)
})
test_that("test_same_crit", {
  expect_equal(fit.mgcca$crit, fit.rgcca$crit)
})
