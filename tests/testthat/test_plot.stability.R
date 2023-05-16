#' # plot.stability
#'''
set.seed(0)
data(Russett)
blocks <- list(
  agriculture = Russett[, seq(3)],
  industry = Russett[, 4:5],
  politic = Russett[, 6:11]
)
fit.sgcca <- rgcca(blocks, sparsity = c(.8, .9, .6))
res <- rgcca_stability(
  fit.sgcca, n_boot = 10, verbose = FALSE, keep = rep(.1, 3)
)

test_that("plot.stability produces the expected plot", {
  skip_if_not(as.logical(Sys.getenv("TEST_SNAPSHOTS")))
  vdiffr::expect_doppelganger(
    "stability weights", plot.stability(res, type = "weights")
  )
})
