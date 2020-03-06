library(rlang)

test_that(
"rgcca_permutation_default", {data("Russett")
    block <- list(
        agriculture = Russett[, seq(3)],
        industry = Russett[, 4:5],
        politic = Russett[, 6:11] )
    expect_is(rgcca_permutation(blocks, n_cores = 1), "permutation")
}
)