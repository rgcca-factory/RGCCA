g  = function(x)  x^2
dg = Deriv::Deriv(g, env = parent.frame())

### Test update_mgcca for matrix blocks
data(Russett)
X_agric = as.matrix(Russett[,c("gini","farm","rent")]);
X_ind = as.matrix(Russett[,c("gnpr","labo")]);
X_polit = as.matrix(Russett[ , c("demostab", "dictator")]);
A = list(X_agric, X_ind, X_polit);
A = scaling(A, scale = T, bias = T, scale_block = T)
C = 1 - diag(3)

test_that("Test that update_mgcca generate a vector a that satisfies the norm
          constraint for 2D blocks and that criterion increases", {
  res_init_mgcca = init_mgcca(A, A, tau = c(1, 1, 1))
  a = res_init_mgcca$a; factors = res_init_mgcca$factors; XtX = res_init_mgcca$XtX
  Y = matrix(0, nrow(X_agric), 3)
  for (j in 1:3) Y[, j] <- A[[j]] %*% a[[j]]
  crit_old = sum(C * g(cov2(Y, bias = T)))
  res_update_mgcca = update_mgcca(A, A, a, factors, XtX, Y, g, dg, C)
  crit = sum(C * g(cov2(res_update_mgcca$Y, bias = T)))
  for (j in 1:3) {
    expect_equal(
      drop(t(res_update_mgcca$a[[j]]) %*% XtX[[j]] %*% res_update_mgcca$a[[j]]),
      1
    )
  }
  expect_true(crit >= crit_old)

  res_init_mgcca = init_mgcca(A, A, tau = c(0.5, 0.5, 0.5))
  a = res_init_mgcca$a; factors = res_init_mgcca$factors; XtX = res_init_mgcca$XtX
  Y = matrix(0, nrow(X_agric), 3)
  for (j in 1:3) Y[, j] <- A[[j]] %*% a[[j]]
  crit_old = sum(C * g(cov2(Y, bias = T)))
  res_update_mgcca = update_mgcca(A, A, a, factors, XtX, Y, g, dg, C)
  crit = sum(C * g(cov2(res_update_mgcca$Y, bias = T)))
  for (j in 1:3) {
    expect_equal(
      drop(t(res_update_mgcca$a[[j]]) %*% XtX[[j]] %*% res_update_mgcca$a[[j]]),
      1
    )
  }
  expect_true(crit >= crit_old)

  res_init_mgcca = init_mgcca(A, A, tau = c(0, 0, 0))
  a = res_init_mgcca$a; factors = res_init_mgcca$factors; XtX = res_init_mgcca$XtX
  Y = matrix(0, nrow(X_agric), 3)
  for (j in 1:3) Y[, j] <- A[[j]] %*% a[[j]]
  crit_old = sum(C * g(cov2(Y, bias = T)))
  res_update_mgcca = update_mgcca(A, A, a, factors, XtX, Y, g, dg, C)
  crit = sum(C * g(cov2(res_update_mgcca$Y, bias = T)))
  for (j in 1:3) {
    expect_equal(
      drop(t(res_update_mgcca$a[[j]]) %*% XtX[[j]] %*% res_update_mgcca$a[[j]]),
      1
    )
  }
  expect_true(crit >= crit_old)
})

### Test update_mgcca for tensor blocks
A = helper.generate_blocks(list(
  c(40, 20, 30), c(40, 35), c(40, 18, 25, 7)
))
A = scaling(A, scale = T, bias = T, scale_block = T)
A_m = lapply(1:3, function(x) matrix(as.vector(A[[x]]), nrow = 40))
C = 1 - diag(3)

test_that("Test that update_mgcca generate a vector a that satisfies the norm
          constraint for tensor blocks, factors are orthogonal and criterion
          increases", {
            res_init_mgcca = init_mgcca(A, A_m, tau = c(1, 1, 1), ranks = c(3, 3, 3))
            a = res_init_mgcca$a; factors = res_init_mgcca$factors; XtX = res_init_mgcca$XtX
            Y = matrix(0, nrow(A[[1]]), 3)
            for (j in 1:3) Y[, j] <- A_m[[j]] %*% a[[j]]
            crit_old = sum(C * g(cov2(Y, bias = T)))
            res_update_mgcca = update_mgcca(A, A_m, a, factors, XtX, Y, g, dg, C, ranks = c(3, 3, 3))
            crit = sum(C * g(cov2(res_update_mgcca$Y, bias = T)))
            for (j in 1:3) {
              expect_equal(
                drop(t(res_update_mgcca$a[[j]]) %*% XtX[[j]] %*% res_update_mgcca$a[[j]]),
                1
              )
              if (length(dim(A[[j]])) > 2) {
                W = Reduce("khatri_rao", rev(res_update_mgcca$factors[[j]]))
                P = t(W) %*% XtX[[j]] %*% W
                diag(P) = 0
                expect_true(max(abs(P)) < 1e-12)
              }
            }
            expect_true(crit >= crit_old)

            res_init_mgcca = init_mgcca(A, A_m, tau = c(0.5, 0.5, 0.5), ranks = c(3, 3, 3))
            a = res_init_mgcca$a; factors = res_init_mgcca$factors; XtX = res_init_mgcca$XtX
            Y = matrix(0, nrow(A[[1]]), 3)
            for (j in 1:3) Y[, j] <- A_m[[j]] %*% a[[j]]
            crit_old = sum(C * g(cov2(Y, bias = T)))
            res_update_mgcca = update_mgcca(A, A_m, a, factors, XtX, Y, g, dg, C, ranks = c(3, 3, 3))
            crit = sum(C * g(cov2(res_update_mgcca$Y, bias = T)))
            for (j in 1:3) {
              expect_equal(
                drop(t(res_update_mgcca$a[[j]]) %*% XtX[[j]] %*% res_update_mgcca$a[[j]]),
                1
              )
              if (length(dim(A[[j]])) > 2) {
                W = Reduce("khatri_rao", rev(res_update_mgcca$factors[[j]]))
                P = t(W) %*% XtX[[j]] %*% W
                diag(P) = 0
                expect_true(max(abs(P)) < 1e-12)
              }
            }
            expect_true(crit >= crit_old)

            # Higher tolerance here as AtA is singular
            res_init_mgcca = init_mgcca(A, A_m, tau = c(0, 0, 0), ranks = c(3, 3, 3))
            a = res_init_mgcca$a; factors = res_init_mgcca$factors; XtX = res_init_mgcca$XtX
            Y = matrix(0, nrow(A[[1]]), 3)
            for (j in 1:3) Y[, j] <- A_m[[j]] %*% a[[j]]
            crit_old = sum(C * g(cov2(Y, bias = T)))
            res_update_mgcca = update_mgcca(A, A_m, a, factors, XtX, Y, g, dg, C, ranks = c(3, 3, 3))
            crit = sum(C * g(cov2(res_update_mgcca$Y, bias = T)))
            for (j in 1:3) {
              expect_equal(
                drop(t(res_update_mgcca$a[[j]]) %*% XtX[[j]] %*% res_update_mgcca$a[[j]]),
                1
              )
              if (length(dim(A[[j]])) > 2) {
                W = Reduce("khatri_rao", rev(res_update_mgcca$factors[[j]]))
                P = t(W) %*% XtX[[j]] %*% W
                diag(P) = 0
                expect_true(max(abs(P)) < 1e-6)
              }
            }
            expect_true(crit >= crit_old)
          })
