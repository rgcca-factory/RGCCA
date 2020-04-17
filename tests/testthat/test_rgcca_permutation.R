data("Russett")
blocks <- list(
    agriculture = Russett[, seq(3)],
    industry = Russett[, 4:5],
    politic = Russett[, 6:11] )

test_that("rgcca_permutation_default", {
        expect_is(rgcca_permutation(blocks, n_cores = 1), "permutation")
    }
)

# test_that("rgcca_permutation_optimal_tau", {
#     expect_is(rgcca_permutation(blocks, tau = "optimal", n_cores = 1), "permutation")
# }
# )