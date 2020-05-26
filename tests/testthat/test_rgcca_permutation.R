data("Russett")
blocks <- list(
    agriculture = Russett[, seq(3)],
    industry = Russett[, 4:5],
    politic = Russett[, 6:11] )

res=rgcca_permutation(blocks, n_cores = 1,nperm = 21)

res=rgcca_permutation(blocks, n_cores = 1,nperm = 10)

test_that("rgcca_permutation_default", {
        expect_is(res, "permutation")
    }
)

# M=matrix(c(0.6,0.6,0.8,0.85,0.7,0.8,0.8,0.9), 2,4)
# 
# res=rgcca_permutation(blocks, n_cores = 1,type="sgcca",superblock=TRUE,perm.par="sparsity")
# res=rgcca_permutation(blocks, n_cores = 1,superblock=TRUE,perm.par="sparsity",perm.value=c(0.8,0.72,0.43,0.5))
# res=rgcca_permutation(blocks, n_cores = 1,type="sgcca",superblock=TRUE,perm.par="sparsity",perm.value=M)
# res=rgcca_permutation(blocks, n_cores = 1,superblock=TRUE,perm.par="tau")
# res=rgcca_permutation(blocks, n_cores = 1,superblock=TRUE,perm.par="tau",perm.value=c(0.8,0.72,0.43,0.5))
# res=rgcca_permutation(blocks, n_cores = 1,superblock=TRUE,perm.par="tau",perm.value=M)
# 
# res=rgcca_permutation(blocks, n_cores = 1,type="sgcca",superblock=FALSE,perm.par="sparsity")
# res=rgcca_permutation(blocks, n_cores = 1,superblock=FALSE,perm.par="sparsity",perm.value=c(0.8,0.72,0.43))
# res=rgcca_permutation(blocks, n_cores = 1,type="sgcca",superblock=FALSE,perm.par="sparsity",perm.value=M[,1:3])
# res=rgcca_permutation(blocks, n_cores = 1,superblock=FALSE,perm.par="tau")
# res=rgcca_permutation(blocks, n_cores = 1,superblock=FALSE,perm.par="tau",perm.value=c(0.8,0.72,0.43))
# res=rgcca_permutation(blocks, n_cores = 1,superblock=FALSE,perm.par="tau",perm.value=M[,1:3])


# test_that("rgcca_permutation_optimal_tau", {
#     expect_is(rgcca_permutation(blocks, tau = "optimal", n_cores = 1), "permutation")
# }
# )