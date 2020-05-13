data("Russett")
blocks <- list(
    agriculture = Russett[, seq(3)],
    industry = Russett[, 4:5],
    politic = Russett[, 6:11] )

# ncomp=1
rgcca_out <- rgcca(blocks, superblock = FALSE,ncomp=1)
resb <- bootstrap_k(rgcca_out, blocks, superblock = FALSE)

test_that("test_bootstrapk_1", {
    expect_is(resb, "list")
    expect_is(resb[[1]], "matrix")
    expect_equal(length(resb), 3)
    expect_true(all(sapply(resb, NCOL) == 1))
    expect_identical(sapply(resb, NROW), sapply(blocks, NCOL))
})

# ncomp=2
rgcca_out_2 <- rgcca(blocks, superblock = FALSE,ncomp=2)
resb_2 <- bootstrap_k(rgcca_out_2, lapply(blocks, scale), superblock = FALSE)

test_that("test_bootstrapk", {
        expect_is(resb_2, "list")
        expect_is(resb_2[[1]], "matrix")
        expect_equal(length(resb_2), 3)
        expect_true(all(sapply(resb_2, NCOL) == 2))
        expect_identical(sapply(resb_2, NROW), sapply(blocks, NCOL))
})
