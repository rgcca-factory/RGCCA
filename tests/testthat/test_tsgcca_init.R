verify_norm_constraints = function(A, A_m, C, sparsity, ranks, init = "svd") {
  res_tsgcca_init = tsgcca_init(A, A_m, C, sparsity = sparsity, ranks = ranks, init = init)
  DIM             = lapply(A, dim)
  pjs             = sapply(DIM, function(x) prod(x[-1]))
  const           = sparsity * sqrt(pjs)
  for (j in 1:length(A)) {
    expect_lte(norm(res_tsgcca_init$a[[j]], type = "2"), 1)
    expect_lte(sum(abs(res_tsgcca_init$a[[j]])), const[j])
  }
}

### Test tsgcca_init for matrix blocks
data(Russett)
X_agric = as.matrix(Russett[,c("gini","farm","rent")]);
X_ind = as.matrix(Russett[,c("gnpr","labo")]);
X_polit = as.matrix(Russett[ , 6:11]);
A = list(X_agric, X_ind, X_polit);
A = scaling(A, scale = T, bias = T, scale_block = T)
C = 1 - diag(length(A))
sparsity = c(0.6, 0.8, 0.8)

test_that("Test that tsgcca_init generate a vector a that satisfies the norm constraints for 2D blocks", {
  verify_norm_constraints(A, A, C, sparsity, ranks = c(1, 1, 1), init = "svd")
  verify_norm_constraints(A, A, C, sparsity, ranks = c(1, 1, 1), init = "random")
})

### Test tsgcca_init for tensor blocks
set.seed(0)
A = helper.generate_blocks(list(
  c(40, 20, 30), c(40, 35), c(40, 18, 25, 7)
))
A = scaling(A, scale = T, bias = T, scale_block = T)
A_m = lapply(1:3, function(x) matrix(as.vector(A[[x]]), nrow = 40))
C = 1 - diag(length(A))
sparsity = c(0.6, 0.8, 0.7)

test_that("Test that tsgcca_init generate a vector a that satisfies the norm
          constraints for tensor blocks", {
            verify_norm_constraints(A, A_m, C, sparsity, ranks = c(1, 1, 1), init = "svd")
            verify_norm_constraints(A, A_m, C, sparsity, ranks = c(1, 1, 1), init = "random")

            verify_norm_constraints(A, A_m, C, sparsity, ranks = c(3, 1, 3), init = "svd")
            verify_norm_constraints(A, A_m, C, sparsity, ranks = c(3, 1, 3), init = "random")
          })
