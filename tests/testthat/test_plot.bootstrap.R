#' # test plot.bootstrap
#-----------------------
set.seed(0)
data("Russett")
blocks <- list(
  agriculture = Russett[, seq(3)],
  industry = Russett[, 4:5],
  politic = Russett[, 6:11]
)
fit.rgcca <- rgcca(blocks, ncomp = 2, method = "rgcca", tau = 1)
fit.boot <- bootstrap(fit.rgcca, n_boot = 20, n_cores = 1, verbose = FALSE)

test_that("plot.bootstrap", {
  vdiffr::expect_doppelganger(
    "Bootstrap weight", plot.bootstrap(fit.boot, type = "weight")
  )

  vdiffr::expect_doppelganger(
    "Bootstrap loadings", plot.bootstrap(
      fit.boot, type = "loadings", show_sign = FALSE
    )
  )
})
