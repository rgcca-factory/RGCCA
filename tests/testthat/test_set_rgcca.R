set.seed(1)

data("Russett")
blocks <- list(
    agriculture = Russett[, seq(3)],
    industry = Russett[, 4:5],
    politic = Russett[, 6:11])

test_that("set_rgcca_equal_to_rgcca", {
    rgcca_out <- rgcca(blocks, response = 2)
    res <- set_rgcca(
        rgcca_out,
        inds =  .Machine$integer.max,
        blocks = blocks,
        response = 2,
        tol = 1E-8
    )
    expect_identical(round(res$Y[[1]], 7), round(rgcca_out$Y[[1]], 7))
})