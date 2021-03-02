# '# Test mgcca
# 
# '''
set.seed(0)

# Generating fake data
A = helper.generate_blocks(list(
  c(40, 20, 30), c(40, 35), c(40, 18, 25, 7)
))
C  = matrix(1, nrow = 3, ncol = 3) - diag(3)

# Run directly mgcca or through calling rgcca. Both should run without errors.
fit.mgcca = mgcca(A, C, ncomp=rep(3, 3, 3), init = "svd", tau = c(1, 1, 1),
                  scheme = "factorial", scale = TRUE, verbose=FALSE,
                  ranks = c(3, 1, 2))
fit.rgcca = rgcca(blocks = A, connection = C, ncomp=rep(3, 3, 3), init = "svd", 
                  tau = c(1, 1, 1), scheme = "factorial", scale = TRUE, 
                  verbose=FALSE, type = "mgcca", ranks = c(3, 1, 2))

# Check number, dimensions and norms of factors
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
test_that("test_norm_factors", {
  # For each component, the sum of square norms of factors of the 1st mode 
  # should equal 1
  expect_equal(sum(diag(crossprod(fit.mgcca$factors[[1]][[1]][, 1:3]))), 1)
  expect_equal(sum(diag(crossprod(fit.mgcca$factors[[1]][[1]][, 4:6]))), 1)
  expect_equal(sum(diag(crossprod(fit.mgcca$factors[[1]][[1]][, 7:9]))), 1)
  expect_equal(sum(diag(crossprod(fit.mgcca$factors[[3]][[1]][, 1:2]))), 1)
  expect_equal(sum(diag(crossprod(fit.mgcca$factors[[3]][[1]][, 3:4]))), 1)
  expect_equal(sum(diag(crossprod(fit.mgcca$factors[[3]][[1]][, 5:6]))), 1)
  
  # For each component, the other factor norms should equal 1
  expect_equal(diag(crossprod(fit.mgcca$factors[[1]][[2]])), rep(1, 3 * 3))
  expect_equal(diag(crossprod(fit.mgcca$factors[[3]][[2]])), rep(1, 3 * 2))
  expect_equal(diag(crossprod(fit.mgcca$factors[[3]][[3]])), rep(1, 3 * 2))
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

# MGCCA should run without errors with random initialization
fit.mgcca = rgcca(blocks = A, connection = C, ncomp=rep(3, 3, 3), 
                  init = "random", tau = c(1, 1, 1), scheme = "factorial", 
                  scale = TRUE, verbose=FALSE, type = "mgcca", 
                  ranks = c(3, 1, 2))

# MGCCA should run without errors with well formatted regularization matrices
regularisation_matrices = list(
  list(
    diag(20), diag(30)
  ), NULL, NULL
)
fit.mgcca = rgcca(blocks = A, connection = C, ncomp=rep(3, 3, 3), 
                  init = "random", tau = c(1, 1, 1), scheme = "factorial", 
                  scale = TRUE, verbose=FALSE, type = "mgcca", 
                  ranks = c(3, 1, 2),
                  regularisation_matrices = regularisation_matrices)

# MGCCA should throw errors when regularization matrices are not well formatted
test_that("test_reg_not_list", {
  expect_error(rgcca(blocks = A, connection = C, 
                     regularisation_matrices = diag(20)), 
               "regularisation_matrices must be NULL or a list of list of matrices", 
               fixed=TRUE)
})
test_that("test_reg_not_list_of_lists", {
  expect_error(rgcca(blocks = A, connection = C, 
                     regularisation_matrices = list(diag(20))), 
               "regularisation_matrices[[1]] must be NULL 
        or a list of matrices", 
               fixed=TRUE)
})
test_that("test_reg_not_square_matrices", {
  expect_error(rgcca(blocks = A, connection = C, 
                     regularisation_matrices = list(list(
    matrix(1, 10, 20)
  ))), 
  "regularisation_matrices matrices must be square matrices", 
  fixed=TRUE)
})
test_that("test_reg_not_enough_matrices", {
  expect_error(rgcca(blocks = A, connection = C, 
                     regularisation_matrices = list(list(
    matrix(1, 20, 20)
  ))), 
  "There should be as many regularisation_matrices 
                          matrices as modes in the block. Mismatch found for
                          block 1.", 
  fixed=TRUE)
})
test_that("test_reg_not_matching_dims", {
  expect_error(rgcca(blocks = A, connection = C, 
                     regularisation_matrices = list(list(
    matrix(1, 20, 20), matrix(-1, 10, 10)
  ))), 
  "regularisation_matrices matrices should match the 
                          mode dimensions. Mismatch found for block 1.", 
  fixed=TRUE)
})

# MGCCA should throw an understandable error if a regularization matrix 
# is singular
test_that("test_reg_not_invertible", {
  expect_error(rgcca(blocks = A, connection = C, 
                     regularisation_matrices = list(list(
    matrix(0, 20, 20), diag(30)
  ))), 
  "Mode 1 regularization matrix for block 1 is singular, please give an invertible matrix.", 
  fixed=TRUE)
})

# MGCCA and RGCCA should perform similarly on matrix blocks
data(Russett)
X_agric = as.matrix(Russett[,c("gini","farm","rent")])
X_ind   = as.matrix(Russett[,c("gnpr","labo")])
X_polit = as.matrix(Russett[ , c("demostab", "dictator")])
A       = list(X_agric, X_ind, X_polit)

fit.rgcca = rgcca(blocks=A, type = "rgcca", connection = 1 - diag(3),
                  scheme = "factorial", tau = "optimal")
fit.mgcca = rgcca(blocks=A, type = "mgcca", connection = 1 - diag(3),
                  scheme = "factorial", tau = "optimal")

test_that("test_highly_correlated_results", {
  for (j in 1:length(A)) {
    expect_true(abs(cor(fit.rgcca$a[[j]], fit.mgcca$a[[j]])) > 1 - 1e-04)
  }
})
