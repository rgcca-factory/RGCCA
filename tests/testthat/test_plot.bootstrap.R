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

test_that("plot.bootstrap raises an error if block is not \"all\" or an integer
          between 1 and the number of blocks", {
  expect_error(
    plot.bootstrap(fit.boot, block = "toto"),
    "block must be equal to \"all\" or an integer between 1 and 3.",
    fixed = TRUE
  )
  expect_error(
    plot.bootstrap(fit.boot, block = 5),
    "block must be equal to \"all\" or an integer between 1 and 3.",
    fixed = TRUE
  )
})

test_that("plot.bootstrap", {
  vdiffr::expect_doppelganger(
    "Bootstrap weight", plot.bootstrap(fit.boot, type = "weight")
  )

  vdiffr::expect_doppelganger(
    "Bootstrap loadings", plot.bootstrap(
      fit.boot, type = "loadings", show_sign = FALSE, block = 1
    )
  )
})
