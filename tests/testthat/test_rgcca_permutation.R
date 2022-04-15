data("Russett")
blocks <- list(
  agriculture = Russett[, seq(3)],
  industry = Russett[, 4:5],
  politic = Russett[, 6:11]
)

res <- rgcca_permutation(blocks, par_type = "tau", n_cores = 1, n_perms = 5)
# res_rgcca=rgcca(blocks)
# res=rgcca_permutation(rgcca_res = res_rgcca, n_cores = 1,n_perms = 21)

res <- rgcca_permutation(blocks, par_type = "tau", n_cores = 1, n_perms = 10)

test_that("rgcca_permutation_default", {
  expect_is(res, "permutation")
})

res_tau <- rgcca_permutation(blocks, par_type = "tau", n_perms = 5, n_cores = 1)
res_sparsity <- rgcca_permutation(blocks, par_type = "sparsity", n_perms = 5, n_cores = 1)

res_tau <- rgcca_permutation(blocks, par_type = "tau", par_value = c(0.7, 0.8, 0.8), n_cores = 1)
res_sparsity <- rgcca_permutation(blocks, par_type = "sparsity", par_value = c(0.8, 0.8, 0.8), n_cores = 1)

test_that("rgcca_sparsity", {
  expect_true(sum(dim(res_sparsity$penalties) == c(10, 3)) == 2)
})
test_that("rgcca_sparsity", {
  expect_true(sum(dim(res_tau$penalties) == c(10, 3)) == 2)
})

M <- matrix(c(0.6, 0.6, 0.8, 0.85, 0.7, 0.8, 0.8, 0.9), 2, 4)
#
# res=rgcca_permutation(blocks, n_cores = 1,method="sgcca",superblock=TRUE,par_type="sparsity")
# res=rgcca_permutation(blocks, n_cores = 1,superblock=TRUE,par_type="sparsity",par_value=c(0.8,0.72,0.43,0.5))
res <- rgcca_permutation(blocks, n_cores = 1, method = "sgcca", superblock = TRUE, par_type = "sparsity", par_value = M)
# res=rgcca_permutation(blocks, n_cores = 1,superblock=TRUE,par_type="tau")
# res=rgcca_permutation(blocks, n_cores = 1,superblock=TRUE,par_type="tau",par_value=c(0.8,0.72,0.43,0.5))
# res=rgcca_permutation(blocks, n_cores = 1,superblock=TRUE,par_type="tau",par_value=M)
#
# res=rgcca_permutation(blocks, n_cores = 1,method="sgcca",superblock=FALSE,par_type="sparsity")
# res=rgcca_permutation(blocks, n_cores = 1,superblock=FALSE,par_type="sparsity",par_value=c(0.8,0.72,0.43))
# res=rgcca_permutation(blocks, n_cores = 1,method="sgcca",superblock=FALSE,par_type="sparsity",par_value=M[,1:3])
# res=rgcca_permutation(blocks, n_cores = 1,superblock=FALSE,par_type="tau")
# res=rgcca_permutation(blocks, n_cores = 1,superblock=FALSE,par_type="tau",par_value=c(0.8,0.72,0.43))
# res=rgcca_permutation(blocks, n_cores = 1,superblock=FALSE,par_type="tau",par_value=M[,1:3])

res_perm <- rgcca_permutation(blocks = blocks, n_cores = 1, par_type = "tau")
test_that("rgcca_crit_perm", {
  expect_true(round(res_perm$crit[length(res_perm$crit)], digits = 2) == 2.42)
})
test_that("rgcca_mat_connection0", {
  expect_true(dim(res_perm$call$connection)[1] == 3)
})
res <- rgcca(res_perm)

res_perm <- rgcca_permutation(blocks = blocks, n_cores = 1, par_type = "tau", superblock = TRUE)
res <- rgcca(res_perm)


# test_that("rgcca_permutation_optimal_tau", {
#     expect_is(rgcca_permutation(blocks, tau = "optimal", n_cores = 1), "permutation")
# }
# )
