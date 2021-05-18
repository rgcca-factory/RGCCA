g  = function(x)  x^2
dg = Deriv::Deriv(g, env = parent.frame())

verify_norm_constraint = function(res, XtX) {
  for (j in 1:length(XtX)) {
    expect_equal(drop(t(res$a[[j]]) %*% XtX[[j]] %*% res$a[[j]]), 1)
  }
}

verify_orthogonality_constraints = function(res, XtX, tol = 1e-12) {
  for (j in 1:length(XtX)) {
    if (length(res$factors[[j]]) > 0) {
      W = Reduce("khatri_rao", rev(res$factors[[j]]))
      P = t(W) %*% XtX[[j]] %*% W
      diag(P) = 0
      expect_true(max(abs(P)) < tol)
    }
  }
}

verify_all = function(A, A_m, tau, ranks, init = "svd", tol = 1e-12) {
  J = length(A)
  C = 1 - diag(J)
  res_ns_mgcca_init = ns_mgcca_init(A, A_m, tau = tau, ranks = ranks, init = init)
  a = res_ns_mgcca_init$a; factors = res_ns_mgcca_init$factors
  weights = res_ns_mgcca_init$weights; XtX = res_ns_mgcca_init$XtX
  Y = matrix(0, nrow(A[[1]]), J)
  for (j in 1:J) Y[, j] <- A_m[[j]] %*% a[[j]]
  crit_old = sum(C * g(cov2(Y, bias = T)))
  res_ns_mgcca_update = ns_mgcca_update(A, A_m, a, factors, weights, XtX, Y, g, dg, C, ranks = ranks)
  crit = sum(C * g(cov2(res_ns_mgcca_update$Y, bias = T)))
  verify_norm_constraint(res_ns_mgcca_update, XtX)
  verify_orthogonality_constraints(res_ns_mgcca_update, XtX, tol = tol)
  expect_true(crit >= crit_old)
}

### Test ns_mgcca_update for matrix blocks
data(Russett)
X_agric = as.matrix(Russett[,c("gini","farm","rent")]);
X_ind = as.matrix(Russett[,c("gnpr","labo")]);
X_polit = as.matrix(Russett[ , c("demostab", "dictator")]);
A = list(X_agric, X_ind, X_polit);
A = scaling(A, scale = T, bias = T, scale_block = T)

test_that("Test that ns_mgcca_update generate a vector a that satisfies the norm
          constraint for 2D blocks and that criterion increases", {
            verify_all(A, A, tau = c(1, 1, 1), ranks = c(1, 1, 1), init = "svd")
            verify_all(A, A, tau = c(1, 1, 1), ranks = c(1, 1, 1), init = "random")

            verify_all(A, A, tau = c(0.5, 0.5, 0.5), ranks = c(1, 1, 1), init = "svd")
            verify_all(A, A, tau = c(0.5, 0.5, 0.5), ranks = c(1, 1, 1), init = "random")

            verify_all(A, A, tau = c(0, 0, 0), ranks = c(1, 1, 1), init = "svd")
            verify_all(A, A, tau = c(0, 0, 0), ranks = c(1, 1, 1), init = "random")
          })

### Test ns_mgcca_update for tensor blocks
set.seed(0)
A = helper.generate_blocks(list(
  c(40, 20, 30), c(40, 35), c(40, 18, 25, 7)
))
A = scaling(A, scale = T, bias = T, scale_block = T)
A_m = lapply(1:3, function(x) matrix(as.vector(A[[x]]), nrow = 40))

test_that("Test that ns_mgcca_update generate a vector a that satisfies the norm
          constraint for tensor blocks, factors are orthogonal and criterion
          increases", {
            verify_all(A, A_m, tau = c(1, 1, 1), ranks = c(3, 3, 3), init = "svd")
            verify_all(A, A_m, tau = c(1, 1, 1), ranks = c(3, 3, 3), init = "random")

            verify_all(A, A_m, tau = c(0.5, 0.5, 0.5), ranks = c(3, 3, 3), init = "svd")
            verify_all(A, A_m, tau = c(0.5, 0.5, 0.5), ranks = c(3, 3, 3), init = "random")

            verify_all(A, A_m, tau = c(0.001, 0.001, 0.001), ranks = c(3, 3, 3), init = "svd", tol = 1e-6)
            verify_all(A, A_m, tau = c(0.001, 0.001, 0.001), ranks = c(3, 3, 3), init = "random", tol = 1e-6)
          })
