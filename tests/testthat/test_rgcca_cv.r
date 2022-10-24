set.seed(1)
data("Russett")
blocks <- list(
  agriculture = Russett[, seq(3)],
  industry = Russett[, 4:5],
  politic = Russett[, 6:8]
)
fit_rgcca <- rgcca(blocks, response = 3)

blocks_classif <- list(
  agriculture = Russett[, 1:3],
  industry = Russett[, 4:5],
  politic = matrix(Russett[, 11], ncol = 1)
)
blocks_classif[["politic"]][blocks_classif[["politic"]][, 1] == 1, ] <- "demo"
blocks_classif[["politic"]][blocks_classif[["politic"]][, 1] == 0, ] <- "ndemo"

test_that("rgcca_cv raises an error if no response is given", {
  expect_error(rgcca_cv(blocks),
    "response is required for rgcca_cv",
    fixed = TRUE
  )
})

test_that("rgcca_cv raises an error if an unknown prediction model is given", {
  expect_error(rgcca_cv(blocks, response = 3, prediction_model = "toto"),
    "unknown model.",
    fixed = TRUE
  )
})

test_that("rgcca_cv raises an error if a regression model is used for a
          classification task", {
  expect_error(rgcca_cv(blocks_classif, response = 3, prediction_model = "lm"),
    "inadequate model.",
    fixed = TRUE
  )
})

test_that("rgcca_cv raises an error if a classification model is used for a
          regression task", {
  expect_error(rgcca_cv(blocks, response = 3, prediction_model = "lda"),
               "inadequate model.",
               fixed = TRUE
  )
})

test_that("rgcca_cv generates a warning if tau is null and block has more
          columns than rows", {
  bad_block <- matrix(rnorm(47 * 127), nrow = 47)
  bad_blocks <- list(bad_block, blocks[[1]], blocks[[2]], blocks[[3]])
  expect_warning(rgcca_cv(bad_blocks, response = 4, prediction_model = "lm"),
    "overfitting risk.",
    fixed = TRUE
  )
  expect_warning(rgcca_cv(blocks, response = 3, prediction_model = "lm"), NA)
})

test_that("rgcca_cv changes par_type to sparsity if a sparse method is given
          with par_type = 'tau'", {
  res <- rgcca_cv(blocks,
    response = 3, par_type = "tau", method = "sgcca",
    par_length = 1, n_run = 1
  )
  expect_equal(res$call$par_type, "sparsity")
  res <- rgcca_cv(blocks,
    response = 3, par_type = "ncomp", method = "sgcca",
    par_length = 1, n_run = 1
  )
  expect_equal(res$call$par_type, "ncomp")
})

test_that("rgcca_cv computes k * n_run scores per parameter value", {
  res <- rgcca_cv(blocks,
    response = 3, par_type = "tau", par_length = 5,
    n_run = 3, validation = "kfold", k = 4
  )
  expect_equal(dim(res$cv), c(5, 4 * 3))
  res <- rgcca_cv(blocks,
    response = 3, par_type = "sparsity", par_length = 2,
    n_run = 5, validation = "kfold", k = 4
  )
  expect_equal(dim(res$cv), c(2, 4 * 5))
  res <- rgcca_cv(blocks,
    response = 3, par_type = "ncomp", par_length = 2,
    n_run = 4, validation = "kfold", k = 7
  )
  expect_equal(dim(res$cv), c(2, 7 * 4))
  res <- rgcca_cv(blocks_classif,
    response = 3, par_type = "ncomp", par_length = 2,
    n_run = 5, validation = "kfold", k = 6, prediction_model = "lda"
  )
  expect_equal(dim(res$cv), c(2, 6 * 5))
  res <- rgcca_cv(fit_rgcca,
    response = 3, par_type = "ncomp", par_length = 1,
    validation = "loo"
  )
  expect_equal(dim(res$cv), c(1, 47))
})
