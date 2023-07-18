#' # summary.rgcca_stability
#'''
set.seed(0)
data(Russett)
blocks <- list(
  agriculture = Russett[, seq(3)],
  industry = Russett[, 4:5],
  politic = Russett[, 6:11]
)
fit.sgcca <- rgcca(blocks, sparsity = c(.8, .9, .6))

test_that("summary.rgcca_stability produces the expected text", {
  skip_if_not(as.logical(Sys.getenv("TEST_SNAPSHOTS")))
  local_edition(3)
  expect_snapshot({
    res <- rgcca_stability(fit.sgcca, n_boot = 10, verbose = FALSE)
    summary(res)
  })
})
