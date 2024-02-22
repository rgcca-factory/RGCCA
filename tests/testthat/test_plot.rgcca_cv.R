#' # test plot.rgcca_cv
#------------------
set.seed(0)
data("Russett")
blocks <- list(
  agriculture = Russett[, seq(3)],
  industry = Russett[, 4:5],
  politic = Russett[, 6:8]
)
res <- rgcca_cv(blocks,
  verbose = FALSE, metric = "MAE",
  response = 3, method = "rgcca", par_type = "tau",
  par_value = c(0, 0.2, 0.3), n_run = 1, n_cores = 1
)
res2 <- suppressWarnings(rgcca_cv(blocks,
  verbose = FALSE,
  response = 3, method = "sgcca", par_type = "ncomp",
  par_length = 3, n_run = 1, n_cores = 1
))

test_that("plot.rgcca_cv produces the expected quantile plot", {
  skip_if_not(as.logical(Sys.getenv("TEST_SNAPSHOTS")))
  vdiffr::expect_doppelganger(
    "CV quantile", plot.rgcca_cv(res, type = "quantile")
  )
})

test_that("plot.rgcca_cv produces the expected sd plot", {
  skip_if_not(as.logical(Sys.getenv("TEST_SNAPSHOTS")))
  vdiffr::expect_doppelganger(
    "CV sd", plot.rgcca_cv(res, type = "sd", display_order = FALSE)
  )
})

test_that("plot.rgcca_cv produces the expected sd plot with many blocks", {
  skip_if_not(as.logical(Sys.getenv("TEST_SNAPSHOTS")))
  vdiffr::expect_doppelganger(
    "CV many blocks", plot.rgcca_cv(res2, type = "sd")
  )
})
