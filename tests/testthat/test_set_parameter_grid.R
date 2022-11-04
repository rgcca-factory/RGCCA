### Generate blocks
data(Russett)
blocks <- list(
  agriculture = Russett[, seq(3)],
  industry = Russett[, 4:5],
  politic = Russett[, 6:11]
)

ncols <- vapply(blocks, NCOL, FUN.VALUE = integer(1))
min_sparsity <- 1 / sqrt(ncols)

test_tau <- function(res, par_length, superblock = FALSE) {
  expect_equal(res$par_type, "tau")
  expect_true(all(res$par_value >= 0))
  expect_true(all(res$par_value <= 1))
  expect_equal(NCOL(res$par_value), length(blocks) + superblock)
  expect_lte(NROW(res$par_value), par_length)
}

test_sparsity <- function(res, par_length, superblock = FALSE) {
  expect_equal(res$par_type, "sparsity")
  for (j in seq_along(blocks)) {
    expect_true(all(res$par_value[, j] >= 1 / sqrt(NCOL(blocks[[j]]))))
  }
  expect_true(all(res$par_value <= 1))
  expect_equal(NCOL(res$par_value), length(blocks) + superblock)
  expect_lte(NROW(res$par_value), par_length)
}

test_ncomp <- function(res, par_length, superblock = FALSE) {
  expect_equal(res$par_type, "ncomp")
  expect_true(all(res$par_value >= 1))
  for (j in seq_along(blocks)) {
    expect_true(all(res$par_value[, j] <= NCOL(blocks[[j]])))
  }
  expect_true(all((res$par_value %% 1) == 0))
  expect_equal(NCOL(res$par_value), length(blocks) + superblock)
  expect_lte(NROW(res$par_value), par_length)
}

### Test set_parameter_grid
test_that("set_parameter_grid returns a valid grid of parameters when par_value
          is NULL", {
  res <- set_parameter_grid("tau", 3, NULL, blocks, response = 3)
  test_tau(res, 3)

  res <- set_parameter_grid("tau", 3, NULL, blocks, response = NULL)
  test_tau(res, 3)

  res <- set_parameter_grid("sparsity", 3, NULL, blocks, response = 3)
  test_sparsity(res, 3)

  res <- set_parameter_grid("ncomp", 3, NULL, blocks, response = 3)
  test_ncomp(res, 3)
})

test_that("set_parameter_grid returns a valid grid of parameters when par_value
          is a valid vector", {
  res <- set_parameter_grid("tau", 3, c(0.5, 1, 0.7), blocks, response = 3)
  test_tau(res, 3)

  res <- set_parameter_grid(
    "tau", 2, c(0.5, 1, 0.7, 1), blocks,
    superblock = TRUE
  )
  test_tau(res, 2, superblock = TRUE)

  res <- set_parameter_grid("sparsity", 3, c(0.9, 1, 0.7), blocks, response = 3)
  test_sparsity(res, 3)

  res <- set_parameter_grid("ncomp", 3, c(1, 2, 2), blocks, response = 3)
  test_ncomp(res, 3)
})

test_that("set_parameter_grid returns a valid grid of parameters when par_value
          is a valid grid", {
  tau <- matrix(pmin(1, pmax(0, rnorm(6))), nrow = 2, ncol = 3)
  res <- set_parameter_grid("tau", 3, tau, blocks, response = 3)
  test_tau(res, 3)

  sparsity <- matrix(c(
    0.6, 0.8, 1,
    0.7, 0.9, 1
  ), nrow = 2, ncol = 3, byrow = TRUE)
  res <- set_parameter_grid("sparsity", 3, sparsity, blocks, response = 3)
  test_sparsity(res, 3)

  ncomp <- matrix(c(1, 2, 2, 3, 1, 3), nrow = 2, ncol = 3, byrow = TRUE)
  res <- set_parameter_grid("ncomp", 3, ncomp, blocks, response = 3)
  test_ncomp(res, 3)
  expect_equal(ncomp, res$par_value)
})

test_that("set_parameter_grid raises errors when par_value is not valid", {
  expect_error(set_parameter_grid("ncomp", 3, print, blocks, response = 3),
    paste0(
      "wrong type of input. par_value must be one of the ",
      "following: NULL, a vector, a matrix or a dataframe."
    ),
    fixed = TRUE
  )
  expect_error(
    set_parameter_grid("ncomp", 3, matrix(1, nrow = 2, ncol = 5), blocks,
      response = 3
    ),
    paste0(
      "wrong shape. If par_value is a matrix or a dataframe,",
      "it must have as many columns as there are blocks (i.e. 3)."
    ),
    fixed = TRUE
  )
  expect_error(
    set_parameter_grid("ncomp", 3, "toto", blocks,
      response = 3
    ),
    paste0("must be numeric"),
    fixed = TRUE
  )
  expect_error(set_parameter_grid("ncomp", 3, 5, blocks, response = 3))
  expect_error(set_parameter_grid("tau", 3, -1, blocks, response = 3))
  expect_error(
    set_parameter_grid("sparsity", 3, matrix(0, nrow = 2, ncol = 3),
      blocks,
      response = 3
    )
  )
})
