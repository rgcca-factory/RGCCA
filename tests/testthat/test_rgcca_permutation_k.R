set.seed(0)
data("Russett")
blocks <- list(
  agriculture = Russett[, seq(3)],
  agriculture2 = Russett[, seq(3)]
)
res_rgcca <- rgcca(blocks)
res_sgcca <- rgcca(blocks, method = "sgcca", sparsity = c(.75, .8))
res_rgcca2 <- rgcca(blocks, ncomp = 2)
res_sgcca2 <- rgcca(blocks, ncomp = 2, method = "sgcca", sparsity = c(.75, .8))
inds <- lapply(blocks, function(x) {
  sample(NROW(x))
})

test_that("rgcca_permutation_k gives the same criterion as rgcca if perm
          is FALSE", {
  res_perm <- rgcca_permutation_k(
    res_rgcca$call, inds, perm = FALSE, par_type = "tau",
    par_value = rep(1, length(blocks))
  )
  expect_equal(res_perm, res_rgcca$crit[length(res_rgcca$crit)])

  res_perm <- rgcca_permutation_k(
    res_rgcca2$call, inds, perm = FALSE, par_type = "tau",
    par_value = rep(1, length(blocks))
  )
  crit <- sum(vapply(
    res_rgcca2$crit, function(x) x[length(x)],
    FUN.VALUE = double(1)
  ))
  expect_equal(res_perm, crit)
})

test_that("rgcca_permutation_k gives a smaller criterion than rgcca if perm
          is TRUE", {
  res_perm <- rgcca_permutation_k(
    res_rgcca$call, inds, perm = TRUE, par_type = "tau",
    par_value = rep(1, length(blocks))
  )
  expect_lt(res_perm, res_rgcca$crit[length(res_rgcca$crit)])

  res_perm <- rgcca_permutation_k(
    res_rgcca2$call, inds, perm = TRUE, par_type = "tau",
    par_value = rep(1, length(blocks))
  )
  crit <- sum(vapply(
    res_rgcca2$crit, function(x) x[length(x)],
    FUN.VALUE = double(1)
  ))
  expect_lt(res_perm, crit)
})

test_that("rgcca_permutation_k gives the same criterion as sgcca if perm
          is FALSE", {
  res_perm <- rgcca_permutation_k(
    res_sgcca$call, inds, perm = FALSE,
    par_type = "sparsity", par_value = c(.75, .8)
  )
  expect_equal(res_perm, res_sgcca$crit[length(res_sgcca$crit)])

  res_perm <- rgcca_permutation_k(
    res_sgcca2$call, inds, perm = FALSE,
    par_type = "sparsity", par_value = c(.75, .8)
  )
  crit <- sum(vapply(
    res_sgcca2$crit, function(x) x[length(x)],
    FUN.VALUE = double(1)
  ))
  expect_equal(res_perm, crit)
})

test_that("rgcca_permutation_k gives a smaller criterion than sgcca if perm
          is TRUE", {
  res_perm <- rgcca_permutation_k(
    res_sgcca$call, inds, perm = TRUE,
    par_type = "sparsity", par_value = c(.75, .8)
  )
  expect_lt(res_perm, res_sgcca$crit[length(res_sgcca$crit)])

  res_perm <- rgcca_permutation_k(
    res_sgcca2$call, inds, perm = TRUE,
    par_type = "sparsity", par_value = c(.75, .8)
  )
  crit <- sum(vapply(
    res_sgcca2$crit, function(x) x[length(x)],
    FUN.VALUE = double(1)
  ))
  expect_lt(res_perm, crit)
})
