data("Russett")
blocks <- list(
    agriculture = Russett[, seq(3)],
    industry = Russett[, 4:5],
    politic = Russett[, 6:11] )

# ncomp=1
rgcca_out <- rgcca(blocks, superblock = FALSE,ncomp=1)
resb <- bootstrap_k(rgcca_out)

test_that("test_bootstrapk_1", {
    expect_is(resb, "list")
    expect_is(resb[[1]][[1]], "matrix")
    expect_is(resb[[2]][[1]], "matrix")
    expect_equal(length(resb), 2)
    expect_true(all(sapply(resb, NCOL) == 1))
    expect_identical(sapply(resb[[1]], NROW), sapply(blocks, NCOL))
    expect_identical(sapply(resb[[2]], NROW), sapply(blocks, NCOL))
})

# ncomp=2
rgcca_out_2 <- rgcca(blocks, superblock = FALSE,ncomp=2)
resb_2 <- bootstrap_k(rgcca_out_2)

test_that("test_bootstrapk", {
        expect_is(resb_2, "list")
        expect_is(resb_2[[1]][[1]], "matrix")
        expect_is(resb_2[[2]][[1]], "matrix")
        expect_equal(length(resb_2), 2)
        expect_true(all(sapply(resb_2[[1]], NCOL) == 2))
        expect_identical(sapply(resb_2[[1]], NROW), sapply(blocks, NCOL))
})

# If one bootstrap sample present at least a single variable with null variance,
# bootstrap_k should return the name of the null variance variables in both the
# two lists it returns.
blocks_3                     = blocks
blocks_3$agriculture$rent    = 0
blocks_3$agriculture$rent[1] = 1
rgcca_out_3                  = rgcca(blocks_3, superblock = FALSE, ncomp = 2)
inds                         = c(2, 2:NROW(blocks_3$agriculture))
resb_3                       = bootstrap_k(rgcca_res = rgcca_out_3, inds = inds)

test_that("test_bootstrapk_missing_var_identification", {
    expect_is(resb_3, "list")
    expect_is(resb_3[[1]], "character")
    expect_is(resb_3[[2]], "character")
    expect_equal(resb_3[[1]], "rent")
    expect_equal(resb_3[[2]], "rent")
})
