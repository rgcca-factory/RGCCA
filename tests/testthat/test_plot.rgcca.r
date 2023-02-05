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

test_that("plot.rgcca produces the expected sample plot", {
  skip_on_cran()
  vdiffr::expect_doppelganger(
    "RGCCA sample", plot.rgcca(
      fit.rgcca, type = "sample", block = seq(2), comp = 1
    )
  )
})

test_that("plot.rgcca produces the expected correlation circle", {
  skip_on_cran()
  vdiffr::expect_doppelganger(
    "RGCCA cor_circle", plot.rgcca(
      fit.rgcca, type = "cor_circle", block = 2,
      comp = seq(2), display_blocks = 2
    )
  )
})

test_that("plot.rgcca produces the expected combined plot with sample plot
          and correlation circle", {
  skip_on_cran()
  vdiffr::expect_doppelganger(
    "RGCCA both", plot.rgcca(fit.rgcca, type = "both", block = 1, comp = seq(2))
  )
})

test_that("plot.rgcca produces the expected combined plot with sample plot
          and correlation circle 2", {
  skip_on_cran()
  vdiffr::expect_doppelganger(
    "RGCCA both 2", plot.rgcca(
      fit.rgcca2, type = "both", block = 4, comp = c(1, 4), show_labels = FALSE
    )
  )
})

test_that("plot.rgcca produces the expected AVE plot", {
  skip_on_cran()
  vdiffr::expect_doppelganger(
    "RGCCA ave", plot.rgcca(fit.rgcca, type = "ave")
  )
})

test_that("plot.rgcca produces the expected weight plot", {
  skip_on_cran()
  vdiffr::expect_doppelganger(
    "RGCCA weight", plot.rgcca(fit.rgcca2, type = "weight", block = 4)
  )
})

test_that("plot.rgcca produces the expected loading plot", {
  skip_on_cran()
  vdiffr::expect_doppelganger(
    "RGCCA loadings", plot.rgcca(fit.rgcca, type = "loadings", block = 1)
  )
})

test_that("plot.rgcca produces the expected loading plot 2", {
  skip_on_cran()
  vdiffr::expect_doppelganger(
    "RGCCA loadings 2", plot.rgcca(fit.rgcca2, type = "loadings")
  )
})
