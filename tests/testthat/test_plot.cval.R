#' # test plot.cval
#------------------
set.seed(0)
data("Russett")
blocks <- list(
  agriculture = Russett[, seq(3)],
  industry = Russett[, 4:5],
  politic = Russett[, 6:8]
)
res <- rgcca_cv(blocks,
  verbose = FALSE,
  response = 3, method = "rgcca", par_type = "tau",
  par_value = c(0, 0.2, 0.3), n_run = 1, n_cores = 1
)
blocks2 <- c(blocks, blocks)
names(blocks2) <- NULL
res2 <- suppressWarnings(rgcca_cv(blocks2,
  verbose = FALSE,
  response = 3, method = "sgcca", par_type = "ncomp",
  par_length = 3, n_run = 1, n_cores = 1
))

test_that("plot.cval", {
  vdiffr::expect_doppelganger(
    "CV quantile", plot.cval(res, type = "quantile")
  )

  vdiffr::expect_doppelganger(
    "CV sd", plot.cval(res, type = "sd", display_order = FALSE)
  )

  vdiffr::expect_doppelganger(
    "CV many blocks", plot.cval(res2, type = "sd")
  )
})
