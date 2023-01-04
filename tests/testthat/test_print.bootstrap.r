#' # print.rgcca_bootstrap
set.seed(0)
data(Russett)
blocks <- list(
  agriculture = Russett[, seq(3)],
  industry = Russett[, 4:5]
)

fit.rgcca <- rgcca(blocks, ncomp = c(1, 2), method = "rgcca", tau = 1)

test_that("print.bootstrap prints the expected string", {
  local_edition(3)
  expect_snapshot({
    res <- rgcca_bootstrap(fit.rgcca, n_boot = 5, n_cores = 1, verbose = FALSE)
    print(res)
  })
})

test_that("print.bootstrap prints the expected string 2", {
  local_edition(3)
  expect_snapshot({
    res <- rgcca_bootstrap(fit.rgcca, n_boot = 2, n_cores = 1, verbose = FALSE)
    print(res, type = "loadings", comp = 2, block = 2)
  })
})
