verify_norm_constraint = function(res, XtX) {
  for (j in 1:length(XtX)) {
    expect_equal(drop(t(res$a[[j]]) %*% XtX[[j]] %*% res$a[[j]]), 1)
  }
}

verify_orthogonality_constraints = function(res, XtX, blocks, tol = 1e-12) {
  for (j in 1:length(blocks)) {
    if (length(dim(blocks[[j]])) > 2) {
      W = Reduce("khatri_rao", rev(res$factors[[j]]))
      P = t(W) %*% XtX[[j]] %*% W
      diag(P) = 0
      expect_true(max(abs(P)) < tol)
    }
  }
}

verify_all = function(A, A_m, tau, ranks, init = "svd", tol = 1e-12) {
  res_init_mgcca = init_mgcca(A, A_m, tau = tau, ranks = ranks, init = init)
  verify_norm_constraint(res_init_mgcca, res_init_mgcca$XtX)
  verify_orthogonality_constraints(res_init_mgcca, res_init_mgcca$XtX, A, tol = tol)
}

### Test init_mgcca for matrix blocks
data(Russett)
X_agric = as.matrix(Russett[,c("gini","farm","rent")]);
X_ind = as.matrix(Russett[,c("gnpr","labo")]);
X_polit = as.matrix(Russett[ , c("demostab", "dictator")]);
A = list(X_agric, X_ind, X_polit);
A = scaling(A, scale = T, bias = T, scale_block = T)

test_that("Test that init_mgcca generate a vector a that satisfies the norm constraint for 2D blocks", {
  verify_all(A, A, tau = c(1, 1, 1), ranks = c(1, 1, 1), init = "svd")
  verify_all(A, A, tau = c(1, 1, 1), ranks = c(1, 1, 1), init = "random")

  verify_all(A, A, tau = c(0.5, 0.5, 0.5), ranks = c(1, 1, 1), init = "svd")
  verify_all(A, A, tau = c(0.5, 0.5, 0.5), ranks = c(1, 1, 1), init = "random")

  verify_all(A, A, tau = c(0, 0, 0), ranks = c(1, 1, 1), init = "svd")
  verify_all(A, A, tau = c(0, 0, 0), ranks = c(1, 1, 1), init = "random")
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
            verify_all(A, A_m, tau = c(1, 1, 1), ranks = c(3, 3, 3), init = "svd")
            verify_all(A, A_m, tau = c(1, 1, 1), ranks = c(3, 3, 3), init = "random")

            verify_all(A, A_m, tau = c(0.5, 0.5, 0.5), ranks = c(3, 3, 3), init = "svd")
            verify_all(A, A_m, tau = c(0.5, 0.5, 0.5), ranks = c(3, 3, 3), init = "random")

            verify_all(A, A_m, tau = c(0, 0, 0), ranks = c(3, 3, 3), init = "svd", tol = 1e-6)
            verify_all(A, A_m, tau = c(0, 0, 0), ranks = c(3, 3, 3), init = "random", tol = 1e-6)
          })
