#' plot.rgcca
#'''
data(Russett)
status <- colnames(Russett)[9:11][apply(Russett[, 9:11], 1, which.max)]
X_agric <- as.matrix(Russett[, c("gini", "farm", "rent")])
X_ind <- as.matrix(Russett[, c("gnpr", "labo")])
X_polit <- as.matrix(Russett[, c("demostab", "dictator")])
blocks <- list(X_agric = X_agric, X_ind = X_ind, X_polit = status)
blocks2 <- list(X_agric = X_agric, X_ind = X_ind, X_polit = X_polit)
fit.rgcca <- rgcca(blocks = blocks, response = 3, tau = 1, ncomp = 2)
fit.rgcca2 <- rgcca(blocks = blocks2, superblock = TRUE, tau = 1, ncomp = 4)

test_that("plot.rgcca produces expected errors", {
  expect_error(
    plot.rgcca(fit.rgcca, type = "biplot", comp = 1),
    "please provide different values for comp",
    fixed = TRUE
  )
  expect_error(
    plot.rgcca(fit.rgcca, type = "cor_circle", block = 3),
    "response components are not orthogonal",
    fixed = TRUE
  )
})

test_that("plot.rgcca produces the expected sample plot", {
  skip_if_not(as.logical(Sys.getenv("TEST_SNAPSHOTS")))
  vdiffr::expect_doppelganger(
    "RGCCA sample", plot.rgcca(
      fit.rgcca, type = "sample", block = seq(2), comp = 1
    )
  )
})

test_that("plot.rgcca produces the expected correlation circle", {
  skip_if_not(as.logical(Sys.getenv("TEST_SNAPSHOTS")))
  vdiffr::expect_doppelganger(
    "RGCCA cor_circle", plot.rgcca(
      fit.rgcca, type = "cor_circle", block = 2,
      comp = seq(2), display_blocks = 2
    )
  )
})

test_that("plot.rgcca produces the expected combined plot with sample plot
          and correlation circle", {
  skip_if_not(as.logical(Sys.getenv("TEST_SNAPSHOTS")))
  vdiffr::expect_doppelganger(
    "RGCCA both", plot.rgcca(fit.rgcca, type = "both", block = 1, comp = seq(2))
  )
})

test_that("plot.rgcca produces the expected combined plot with sample plot
          and correlation circle 2", {
  skip_if_not(as.logical(Sys.getenv("TEST_SNAPSHOTS")))
  vdiffr::expect_doppelganger(
    "RGCCA both 2", plot.rgcca(
      fit.rgcca2, type = "both", block = 4,
      comp = c(1, 4), show_var_names = FALSE
    )
  )
})

test_that("plot.rgcca produces the expected AVE plot", {
  skip_if_not(as.logical(Sys.getenv("TEST_SNAPSHOTS")))
  vdiffr::expect_doppelganger(
    "RGCCA ave", plot.rgcca(fit.rgcca, type = "ave")
  )
})

test_that("plot.rgcca produces the expected weight plot", {
  skip_if_not(as.logical(Sys.getenv("TEST_SNAPSHOTS")))
  vdiffr::expect_doppelganger(
    "RGCCA weight", plot.rgcca(fit.rgcca2, type = "weight", block = 4)
  )
})

test_that("plot.rgcca produces the expected loading plot", {
  skip_if_not(as.logical(Sys.getenv("TEST_SNAPSHOTS")))
  vdiffr::expect_doppelganger(
    "RGCCA loadings", plot.rgcca(fit.rgcca, type = "loadings", block = 1)
  )
})

test_that("plot.rgcca produces the expected loading plot 2", {
  skip_if_not(as.logical(Sys.getenv("TEST_SNAPSHOTS")))
  vdiffr::expect_doppelganger(
    "RGCCA loadings 2", plot.rgcca(fit.rgcca2, type = "loadings")
  )
})

test_that("plot.rgcca produces the expected biplot", {
  skip_if_not(as.logical(Sys.getenv("TEST_SNAPSHOTS")))
  vdiffr::expect_doppelganger(
    "RGCCA biplot", plot.rgcca(fit.rgcca, type = "biplot")
  )
})

test_that("plot.rgcca produces the expected biplot 2", {
  skip_if_not(as.logical(Sys.getenv("TEST_SNAPSHOTS")))
  vdiffr::expect_doppelganger(
    "RGCCA biplot 2", plot.rgcca(
      fit.rgcca2, type = "biplot", show_arrows = FALSE
    )
  )
})

test_that("plot.rgcca produces the expected biplot 3", {
  skip_if_not(as.logical(Sys.getenv("TEST_SNAPSHOTS")))
  vdiffr::expect_doppelganger(
    "RGCCA biplot 3", plot.rgcca(
      fit.rgcca2, type = "biplot",
      response = Russett[, 7], show_sample_names = FALSE
    )
  )
})

test_that("plot.rgcca produces the expected plots with selection", {
  skip_if_not(as.logical(Sys.getenv("TEST_SNAPSHOTS")))
  vdiffr::expect_doppelganger(
    "RGCCA biplot sample_select", plot.rgcca(
      fit.rgcca2, type = "biplot",
      response = Russett[, 7], sample_select = list(name = c("France", "Japan"))
    )
  )
  vdiffr::expect_doppelganger(
    "RGCCA cor_circle var_select", plot.rgcca(
      fit.rgcca2, type = "cor_circle",
      var_select = list(value = .5)
    )
  )
  vdiffr::expect_doppelganger(
    "RGCCA weights var_select", plot.rgcca(
      fit.rgcca, type = "weights",
      var_select = list(number = 1)
    )
  )
  vdiffr::expect_doppelganger(
    "RGCCA loadings var_select", plot.rgcca(
      fit.rgcca, type = "loadings",
      var_select = list(number = 2), display_order = FALSE
    )
  )
})

test_that("plot.rgcca produces expected errors with selection", {
  expect_error(
    plot.rgcca(fit.rgcca, type = "biplot", sample_select = list(toto = 1)),
    paste0(
      "sample_select must be NULL or a named list. Possible names are ",
      "'name', 'value', and 'number'."
    ),
    fixed = TRUE
  )
  expect_error(
    plot.rgcca(fit.rgcca2, type = "weight", var_select = list(name = "toto")),
    paste0(
      "Wrong variable name. The names in var_select$name do not all correspond",
      " to existing variable names."
    ),
    fixed = TRUE
  )
  expect_error(
    plot.rgcca(fit.rgcca2, type = "cor_circle", var_select = list(value = 5)),
    "var_select$value must be lower than or equal to 0.891150801400273.",
    fixed = TRUE
  )
  expect_error(
    plot.rgcca(fit.rgcca2, type = "loadings", var_select = list(number = 0)),
    "var_select$number must be higher than or equal to 1",
    fixed = TRUE
  )
})
