#' # test plot.rgcca_bootstrap
#-----------------------
set.seed(0)
data("Russett")
blocks <- list(
  agriculture = Russett[, seq(3)],
  industry = Russett[, 4:5],
  politic = Russett[, 6:11]
)
fit.rgcca <- rgcca(blocks, ncomp = 2, method = "rgcca", tau = 1)
fit.boot <- rgcca_bootstrap(
  fit.rgcca, n_boot = 20, n_cores = 1, verbose = FALSE
)

test_that("plot.rgcca_bootstrap raises an error if block is not an integer
          between 1 and the number of raw blocks", {
  expect_error(
    plot.rgcca_bootstrap(fit.boot, block = "toto"),
    "block must be numeric.",
    fixed = TRUE
  )
  expect_error(
    plot.rgcca_bootstrap(fit.boot, block = 5),
    "block must be lower than the number of blocks, i.e. 3.",
    fixed = TRUE
  )
})

test_that("plot.rgcca_bootstrap produces the expected weight plot", {
  skip_if_not(as.logical(Sys.getenv("TEST_SNAPSHOTS")))
  vdiffr::expect_doppelganger(
    "Bootstrap weight", plot.rgcca_bootstrap(fit.boot, type = "weight")
  )
})

test_that("plot.rgcca_bootstrap produces the expected loading plot", {
  skip_if_not(as.logical(Sys.getenv("TEST_SNAPSHOTS")))
  vdiffr::expect_doppelganger(
    "Bootstrap loadings", plot.rgcca_bootstrap(
      fit.boot, type = "loadings", show_stars = FALSE, block = 1
    )
  )
})
