#' # print.bootstrap
set.seed(0)
data(Russett)
blocks <- list(
  agriculture = Russett[, seq(3)],
  industry = Russett[, 4:5]
)

fit.rgcca <- rgcca(blocks, ncomp = c(1, 2), method = "rgcca", tau = 1)

test_that("print.bootstrap", {
  local_edition(3)
  expect_snapshot({
    res <- bootstrap(fit.rgcca, n_boot = 5, n_cores = 1, verbose = FALSE)
    print(res)
  })

  expect_snapshot({
    res <- bootstrap(fit.rgcca, n_boot = 2, n_cores = 1, verbose = FALSE)
    print(res, type = "loadings")
  })
})