### Test init_mgcca for matrix blocks
data(Russett)
X_agric = as.matrix(Russett[,c("gini","farm","rent")]);
X_ind = as.matrix(Russett[,c("gnpr","labo")]);
X_polit = as.matrix(Russett[ , c("demostab", "dictator")]);
A = list(X_agric, X_ind, X_polit);
A = scaling(A, scale = T, bias = T, scale_block = T)

test_that("Test that init_mgcca generate a vector a that satisfies the norm constraint for 2D blocks", {
  res_init_mgcca = init_mgcca(A, A, tau = c(1, 1, 1))
  for (j in 1:3) {
    expect_equal(
      drop(t(res_init_mgcca$a[[j]]) %*% res_init_mgcca$XtX[[j]] %*% res_init_mgcca$a[[j]]),
      1
    )
  }

  res_init_mgcca = init_mgcca(A, A, tau = c(0.5, 0.5, 0.5))
  for (j in 1:3) {
    expect_equal(
      drop(t(res_init_mgcca$a[[j]]) %*% res_init_mgcca$XtX[[j]] %*% res_init_mgcca$a[[j]]),
      1
    )
  }

  res_init_mgcca = init_mgcca(A, A, tau = c(0, 0, 0))
  for (j in 1:3) {
    expect_equal(
      drop(t(res_init_mgcca$a[[j]]) %*% res_init_mgcca$XtX[[j]] %*% res_init_mgcca$a[[j]]),
      1
    )
  }
})

### Test init_mgcca for tensor blocks
set.seed(0)
A = helper.generate_blocks(list(
  c(40, 20, 30), c(40, 35), c(40, 18, 25, 7)
))
A = scaling(A, scale = T, bias = T, scale_block = T)
A_m = lapply(1:3, function(x) matrix(as.vector(A[[x]]), nrow = 40))

test_that("Test that init_mgcca generate a vector a that satisfies the norm
          constraint for tensor blocks and factors are orthogonal", {
  res_init_mgcca = init_mgcca(A, A_m, tau = c(1, 1, 1), ranks = c(3, 3, 3))
  for (j in 1:3) {
    expect_equal(
      drop(t(res_init_mgcca$a[[j]]) %*% res_init_mgcca$XtX[[j]] %*% res_init_mgcca$a[[j]]),
      1
    )
    if (length(dim(A[[j]])) > 2) {
      W = Reduce("khatri_rao", rev(res_init_mgcca$factors[[j]]))
      P = t(W) %*% res_init_mgcca$XtX[[j]] %*% W
      diag(P) = 0
      expect_true(max(abs(P)) < 1e-12)
    }
  }

  res_init_mgcca = init_mgcca(A, A_m, tau = c(0.5, 0.5, 0.5), ranks = c(3, 3, 3))
  for (j in 1:3) {
    expect_equal(
      drop(t(res_init_mgcca$a[[j]]) %*% res_init_mgcca$XtX[[j]] %*% res_init_mgcca$a[[j]]),
      1
    )
    if (length(dim(A[[j]])) > 2) {
      W = Reduce("khatri_rao", rev(res_init_mgcca$factors[[j]]))
      P = t(W) %*% res_init_mgcca$XtX[[j]] %*% W
      diag(P) = 0
      expect_true(max(abs(P)) < 1e-12)
    }
  }

  # Higher tolerance here as AtA is singular
  res_init_mgcca = init_mgcca(A, A_m, tau = c(0, 0, 0), ranks = c(3, 3, 3))
  for (j in 1:3) {
    expect_equal(
      drop(t(res_init_mgcca$a[[j]]) %*% res_init_mgcca$XtX[[j]] %*% res_init_mgcca$a[[j]]),
      1
    )
    if (length(dim(A[[j]])) > 2) {
      W = Reduce("khatri_rao", rev(res_init_mgcca$factors[[j]]))
      P = t(W) %*% res_init_mgcca$XtX[[j]] %*% W
      diag(P) = 0
      expect_true(max(abs(P)) < 1e-6)
    }
  }
})
