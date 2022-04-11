### Test set_rgcca
set.seed(1)

data("Russett")
blocks <- list(
  agriculture = Russett[, seq(3)],
  industry = Russett[, 4:5], politic = Russett[, 6:11]
)

test_that("set_rgcca gives the same results as rgcca", {
  fit_rgcca <- rgcca(blocks, superblock = TRUE)
  res <- set_rgcca(fit_rgcca)
  expect_identical(fit_rgcca, res)

  fit_rgcca <- rgcca(blocks, response = 2, scale = FALSE, scale_block = FALSE)
  res <- set_rgcca(fit_rgcca)
  expect_identical(fit_rgcca, res)

  inds <- sample(seq(nrow(blocks[[1]])), 10, replace = FALSE)
  fit_rgcca_minus_inds <- rgcca(
    lapply(blocks, function(x) x[-inds, ]),
    response = 2,
    scale = FALSE, scale_block = "lambda1"
  )
  res <- set_rgcca(fit_rgcca, inds = inds, scale_block = "lambda1")
  expect_identical(fit_rgcca_minus_inds, res)
})
