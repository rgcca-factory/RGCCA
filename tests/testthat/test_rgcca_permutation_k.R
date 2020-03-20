data("Russett")
blocks <- list(
    agriculture = Russett[, seq(3)],
    industry = Russett[, 4:5],
    politic = Russett[, 6:11] )

par = expand.grid(rep(list(seq(2)), length(blocks)))

test_structure <- function(res){
    expect_is(res, "numeric")
    expect_equal(length(res), 8)
}

test_that(
    "rgcca_permutationk_default", {
        test_structure(rgcca_permutation_k(blocks, n_cores = 1))
        expect_identical(
            sapply(
                seq(NROW(par)), 
                function(i) 
                    sum(
                        sapply(
                            rgcca(
                                blocks, 
                                tol = 1e-3, 
                                tau = rep(1, 3), 
                                ncomp = par[i, ]
                            )$crit,
                            sum))),
            rgcca_permutation_k(blocks, perm = FALSE, n_cores = 1)
        )
    }
)

test_that(
    "rgcca_permutationk_optimal_tau", {
        test_structure(rgcca_permutation_k(blocks, tau = "optimal", n_cores = 1))
        expect_identical(
            sapply(
                seq(NROW(par)), 
                function(i) 
                    sum(
                        sapply(
                            rgcca(
                                blocks, 
                                tol = 1e-3, 
                                tau = 0, 
                                ncomp = par[i, ]
                            )$crit,
                            sum))),
            rgcca_permutation_k(blocks, tau = 0, perm = FALSE, n_cores = 1)
        )
    }
)