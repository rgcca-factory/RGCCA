set.seed(0)
data("Russett")
blocks <- list(
  agriculture = Russett[, seq(3)],
  agriculture2 = Russett[, seq(3)]
)

test_that("rgcca_permutation_k gives the same criterion as rgcca if perm
          is FALSE", {
  res_perm <- rgcca_permutation_k(blocks, perm = FALSE)
  res_rgcca <- rgcca(blocks)
  expect_equal(res_perm, res_rgcca$crit[length(res_rgcca$crit)])

  res_perm <- rgcca_permutation_k(blocks, perm = FALSE, ncomp = 2)
  res_rgcca <- rgcca(blocks, ncomp = 2)
  crit <- sum(vapply(
    res_rgcca$crit, function(x) x[length(x)],
    FUN.VALUE = double(1)
  ))
  expect_equal(res_perm, crit)
})

test_that("rgcca_permutation_k gives a smaller criterion than rgcca if perm
          is TRUE", {
  res_perm <- rgcca_permutation_k(blocks, perm = TRUE)
  res_rgcca <- rgcca(blocks)
  expect_lt(res_perm, res_rgcca$crit[length(res_rgcca$crit)])

  res_perm <- rgcca_permutation_k(blocks, perm = TRUE, ncomp = 2)
  res_rgcca <- rgcca(blocks, ncomp = 2)
  crit <- sum(vapply(
    res_rgcca$crit, function(x) x[length(x)],
    FUN.VALUE = double(1)
  ))
  expect_lt(res_perm, crit)
})

test_that("rgcca_permutation_k gives the same criterion as sgcca if perm
          is FALSE", {
  res_perm <- rgcca_permutation_k(blocks,
    perm = FALSE, method = "sgcca",
    par_value = c(.75, .8), par_type = "sparsity"
  )
  res_sgcca <- rgcca(blocks, method = "sgcca", sparsity = c(.75, .8))
  expect_equal(res_perm, res_sgcca$crit[length(res_sgcca$crit)])

  res_perm <- rgcca_permutation_k(blocks,
    perm = FALSE, method = "sgcca",
    par_value = c(.75, .8), par_type = "sparsity",
    ncomp = 2
  )
  res_sgcca <- rgcca(blocks, ncomp = 2, method = "sgcca", sparsity = c(.75, .8))
  crit <- sum(vapply(
    res_sgcca$crit, function(x) x[length(x)],
    FUN.VALUE = double(1)
  ))
  expect_equal(res_perm, crit)
})

test_that("rgcca_permutation_k gives a smaller criterion than sgcca if perm
          is TRUE", {
  res_perm <- rgcca_permutation_k(blocks,
    perm = TRUE, method = "sgcca",
    par_value = c(.75, .8), par_type = "sparsity"
  )
  res_sgcca <- rgcca(blocks, method = "sgcca", sparsity = c(.75, .8))
  expect_lt(res_perm, res_sgcca$crit[length(res_sgcca$crit)])

  res_perm <- rgcca_permutation_k(blocks,
    perm = TRUE, method = "sgcca",
    par_value = c(.75, .8), par_type = "sparsity",
    ncomp = 2
  )
  res_sgcca <- rgcca(blocks, ncomp = 2, method = "sgcca", sparsity = c(.75, .8))
  crit <- sum(vapply(
    res_sgcca$crit, function(x) x[length(x)],
    FUN.VALUE = double(1)
  ))
  expect_lt(res_perm, crit)
})
