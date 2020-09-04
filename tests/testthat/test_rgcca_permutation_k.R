data("Russett")
blocks <- list(
    agriculture = Russett[, seq(3)],
    industry = Russett[, 4:5],
    politic = Russett[, 6:11] )

par = expand.grid(rep(list(seq(2)), length(blocks)))
rgcca(
    blocks, 
    tol = 1e-3, 
    tau = rep(1, 3), 
    ncomp = par[1, ]
    )

res_perm=rgcca_permutation_k(blocks,tol=1e-3,tau=rep(1,3),ncomp=par[1,])
res=rgcca_permutation_k(blocks,tol=1e-3,tau=rep(1,3),ncomp=par[1,],perm=FALSE)


test_that(
    "rgcca_permutationk_without_perm", {
        expect_equal(round(res,digits=6),1.434355)    
    }
)
test_that(
    "rgcca_permutationk_default", {
        expect_true(res_perm<1.434355)    
    }
)
# test_that(
#     "rgcca_permutationk_optimal_tau", {
#         test_structure(rgcca_permutation_k(blocks, tau = "optimal", superblock = FALSE, n_cores = 1))
#         expect_identical(
#             sapply(
#                 seq(NROW(par)), 
#                 function(i) 
#                     sum(
#                         sapply(
#                             rgcca(
#                                 blocks, 
#                                 tol = 1e-3, 
#                                 tau = 0, 
#                                 ncomp = par[i, ]
#                             )$crit,
#                             sum))),
#             rgcca_permutation_k(blocks, tau = 0, perm = FALSE, n_cores = 1)
#         )
#     }
# )