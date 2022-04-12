data("Russett")
blocks <- list(
  agriculture = Russett[, seq(3)],
  industry = Russett[, 4:5],
  politic = Russett[, 6:11]
)
fit_rgcca <- rgcca(blocks, superblock = TRUE)

test_that("rgcca_permutation raises an error if only one block is given", {
  expect_error(rgcca_permutation(list(blocks[[1]])),
    "wrong number of blocks.",
    fixed = TRUE
  )
})

test_that("rgcca_permutation changes par_type to sparsity if a sparse method is
          given with par_type = 'tau'", {
  res <- rgcca_permutation(blocks,
    response = 3, par_type = "tau",
    method = "sgcca", par_length = 1, n_perms = 1
  )
  expect_equal(res$call$par_type, "sparsity")
  res <- rgcca_permutation(blocks,
    response = 3, par_type = "ncomp",
    method = "sgcca", par_length = 1, n_perms = 1
  )
  expect_equal(res$call$par_type, "ncomp")
})

test_that("rgcca_permutation computes n_perms permuted scores and one non
          permuted score per parameter value", {
  res <- rgcca_permutation(blocks,
    par_type = "tau", par_length = 5,
    n_perms = 3
  )
  expect_equal(dim(res$permcrit), c(5, 3))
  expect_equal(length(res$crit), 5)
  res <- rgcca_permutation(blocks,
    par_type = "sparsity", par_length = 7,
    n_perms = 5
  )
  expect_equal(dim(res$permcrit), c(7, 5))
  expect_equal(length(res$crit), 7)
  res <- rgcca_permutation(blocks,
    par_type = "ncomp", par_length = 2,
    n_perms = 4
  )
  expect_equal(dim(res$permcrit), c(2, 4))
  expect_equal(length(res$crit), 2)
  res <- rgcca_permutation(fit_rgcca,
    par_value = c(0.5, 1, 1, 1),
    par_length = 1, n_perms = 2
  )
  expect_equal(dim(res$permcrit), c(1, 2))
  expect_equal(length(res$crit), 1)
})

test_that("rgcca imports the parameters from a permutation object", {
  res <- rgcca_permutation(blocks,
    par_type = "sparsity", par_length = 5,
    n_perms = 3
  )
  fit_rgcca <- rgcca(res)
  expect_identical(unname(res$bestpenalties), fit_rgcca$call$sparsity)
})
