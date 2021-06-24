set.seed(0)

g  = function(x)  x^2
dg = Deriv::Deriv(g, env = parent.frame())

verify_norm_constraints = function(A, A_m, C, sparsity, ranks, init = "svd", scheme = "factorial") {
  J          = length(A)
  res_init   = tsgcca_init(A, A_m, C, sparsity = sparsity, ranks = ranks, init = init, scheme = scheme)
  a          = res_init$a; factors = res_init$factors; weights = res_init$weights
  DIM        = lapply(A, dim)
  pjs        = sapply(DIM, function(x) prod(x[-1]))
  const      = sparsity * sqrt(pjs)
  Y          = matrix(0, nrow(A[[1]]), J)
  for (j in 1:J) Y[, j] <- A_m[[j]] %*% a[[j]]
  crit_old   = sum(C * g(cov2(Y, bias = T)))
  res_update = tsgcca_update(A, A_m, a, factors, weights, Y, g, dg, C, ranks = ranks, sparsity = sparsity)
  crit       = sum(C * g(cov2(res_update$Y, bias = T)))
  for (j in 1:length(A)) {
    expect_lte(norm(res_update$a[[j]], type = "2"), 1)
    expect_lte(sum(abs(res_update$a[[j]])), const[j])
  }
  expect_gte(crit, crit_old)
}

### Test tsgcca_update for matrix blocks
data(Russett)
X_agric = as.matrix(Russett[,c("gini","farm","rent")]);
X_ind = as.matrix(Russett[,c("gnpr","labo")]);
X_polit = as.matrix(Russett[ , c("demostab", "dictator")]);
A = list(X_agric, X_ind, X_polit);
A = scaling(A, scale = T, bias = T, scale_block = T)
C = 1 - diag(length(A))
sparsity = c(0.6, 0.8, 0.8)

test_that("Test that tsgcca_update generate a vector a that satisfies the norm
          constraint for 2D blocks and that criterion increases", {
            verify_norm_constraints(A, A, C, sparsity = sparsity, ranks = c(1, 1, 1), init = "svd")
            verify_norm_constraints(A, A, C, sparsity = sparsity, ranks = c(1, 1, 1), init = "random")
          })

### Test tsgcca_update for tensor blocks
A = helper.generate_blocks(list(
  c(40, 20, 30), c(40, 35), c(40, 18, 25, 7)
))
A = scaling(A, scale = T, bias = T, scale_block = T)
A_m = lapply(1:3, function(x) matrix(as.vector(A[[x]]), nrow = 40))
C = 1 - diag(length(A))
sparsity = c(0.6, 0.8, 0.7)

test_that("Test that tsgcca_update generate a vector a that satisfies the norm
          constraint for tensor blocks and criterion increases", {
            verify_norm_constraints(A, A_m, C, sparsity, ranks = c(1, 1, 1), init = "svd")
            verify_norm_constraints(A, A_m, C, sparsity, ranks = c(1, 1, 1), init = "random")

            verify_norm_constraints(A, A_m, C, sparsity, ranks = c(3, 1, 3), init = "svd")
            verify_norm_constraints(A, A_m, C, sparsity, ranks = c(3, 1, 3), init = "random")
          })
