data("Russett")
blocks <- list(
  agriculture = Russett[, seq(3)],
  industry = Russett[, 4:5],
  politic = Russett[, 6:11]
)
p <- vapply(blocks, NCOL, FUN.VALUE = 1L)

### Case ncomp = 1
rgcca_out <- rgcca(blocks, superblock = FALSE, ncomp = 1)
resb <- rgcca_bootstrap_k(rgcca_out)

test_that("test_rgcca_bootstrap_k_1", {
  expect_is(resb, "list")
  expect_is(resb[[1]][[1]], "matrix")
  expect_is(resb[[2]][[1]], "matrix")
  expect_equal(length(resb), 2)
  expect_true(all(vapply(resb, NCOL, FUN.VALUE = 1L) == 1))
  expect_identical(vapply(resb[[1]], NROW, FUN.VALUE = 1L), p)
  expect_identical(vapply(resb[[2]], NROW, FUN.VALUE = 1L), p)
})

### Case ncomp = 2
rgcca_out_2 <- rgcca(blocks, superblock = FALSE, ncomp = 2)
resb_2 <- rgcca_bootstrap_k(rgcca_out_2)

test_that("test_rgcca_bootstrap_k_2", {
  expect_is(resb_2, "list")
  expect_is(resb_2[[1]][[1]], "matrix")
  expect_is(resb_2[[2]][[1]], "matrix")
  expect_equal(length(resb_2), 2)
  expect_true(all(vapply(resb_2[[1]], NCOL, FUN.VALUE = 1L) == 2))
  expect_identical(vapply(resb_2[[1]], NROW, FUN.VALUE = 1L), p)
})

# If one bootstrap sample presents at least a single variable with null
# variance, rgcca_bootstrap_k should still return results
blocks_3 <- blocks
blocks_3$agriculture$rent <- 0
blocks_3$agriculture$rent[1] <- 1
rgcca_out_3 <- rgcca(blocks_3, superblock = FALSE, ncomp = 2)
inds <- c(2, 2:NROW(blocks_3$agriculture))
resb_3 <- rgcca_bootstrap_k(rgcca_res = rgcca_out_3, inds = inds)

test_that("test_rgcca_bootstrap_k_3", {
  expect_is(resb_3, "list")
  expect_is(resb_3[[1]][[1]], "matrix")
  expect_is(resb_3[[2]][[1]], "matrix")
  expect_equal(length(resb_3), 2)
  expect_true(all(vapply(resb_3[[1]], NCOL, FUN.VALUE = 1L) == 2))
  expect_identical(vapply(resb_3[[1]], NROW, FUN.VALUE = 1L), p)
})
