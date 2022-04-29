data(Russett)
blocks <- list(
  agriculture = Russett[, seq(3)],
  industry = Russett[, 4:5],
  politic = Russett[, 6:11]
)

test_that("print_comp returns the AVE of the specified component", {
  fit.rgcca <- rgcca(blocks)
  res <- print_comp(fit.rgcca, n = 1, i = 1, outer = FALSE)
  expect_equal(res, "Comp. 1 (69.9%)")

  fit.rgcca <- rgcca(blocks, ncomp = 2, sparsity = c(1, 1, 0.5))
  res <- print_comp(fit.rgcca, n = 2, i = 3, outer = FALSE)
  expect_equal(res, "Comp. 2 (3 variables, 18.2%)")
})

test_that("print_comp returns the outer AVE", {
  fit.rgcca <- rgcca(blocks)
  res <- print_comp(fit.rgcca, outer = TRUE)
  expect_equal(res, "First outer AVE: 60.2%")

  fit.rgcca <- rgcca(blocks, ncomp = 4, superblock = TRUE)
  res <- print_comp(fit.rgcca, outer = TRUE)
  expect_equal(res, "First corrected outer AVE:  50.6% & 12.9%")
})
