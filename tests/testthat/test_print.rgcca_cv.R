#' # print.rgcca_cv
#'''
set.seed(0)
data(Russett)
blocks <- list(
  agriculture = Russett[, seq(3)],
  industry = Russett[, 4:5],
  politic = Russett[, 6:8]
)
blocks_classif <- list(
  agriculture = Russett[, seq(3)],
  industry = Russett[, 4:5],
  politic = as.factor(Russett[, 9])
)

test_that("print.rgcca_cv prints the expected text", {
  skip_if_not(as.logical(Sys.getenv("TEST_SNAPSHOTS")))
  local_edition(3)
  expect_snapshot({
    res <- rgcca_cv(blocks,
      response = 3, method = "rgcca", par_type = "tau",
      par_value = c(0, 0.2, 0.3), n_run = 1, n_cores = 1,
      par_length = 2, verbose = FALSE
    )
    print(res)
  })
})

test_that("print.rgcca_cv prints the expected text 2", {
  skip_if_not(as.logical(Sys.getenv("TEST_SNAPSHOTS")))
  local_edition(3)
  expect_snapshot({
    res <- rgcca_cv(blocks_classif,
                    response = 3, method = "rgcca", par_type = "tau",
                    par_value = c(0, 0.2, 0.3), n_run = 1, n_cores = 1,
                    par_length = 2, verbose = FALSE, prediction_model = "lda"
    )
    print(res)
  })
})
