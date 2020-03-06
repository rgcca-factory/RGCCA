data("Russett")
block <- list(
    agriculture = Russett[, seq(3)],
    industry = Russett[, 4:5],
    politic = Russett[, 6:11] )

rgcca_out <- rgcca(blocks)
bootstrap_k(rgcca_out)
resb <- bootstrap_k(rgcca_out, lapply(blocks, scale), superblock = FALSE)

test_that("test_bootstrapk", {
        expect_is(resb, "list")
        expect_is(resb[[1]], "matrix")
        expect_equal(length(resb), 3)
        expect_true(all(sapply(resb, NCOL) == 2))
        expect_identical(sapply(resb, NROW), sapply(blocks, NCOL))
})
