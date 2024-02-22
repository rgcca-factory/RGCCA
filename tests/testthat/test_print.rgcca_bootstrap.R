#' # print.rgcca_bootstrap
set.seed(0)
data(Russett)
blocks <- list(
  agriculture = Russett[, seq(3)],
  industry = Russett[, 4:5]
)

fit_rgcca <- rgcca(blocks, ncomp = c(1, 2), method = "rgcca", tau = 1)
res <- rgcca_bootstrap(fit_rgcca, n_boot = 5, n_cores = 1, verbose = FALSE)

test_that("print.rgcca_bootstrap prints the expected string", {
  skip_if_not(as.logical(Sys.getenv("TEST_SNAPSHOTS")))
  local_edition(3)
  expect_snapshot(print(res))
})
