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

test_that("rgcca_permutation_k gives the same criterion as rgcca if perm
          is FALSE", {
  res_perm <- rgcca_permutation_k(
    res_rgcca$call, perm = FALSE, par_type = "tau",
    par_value = rep(1, length(blocks))
  )$crit
  expect_equal(res_perm, res_rgcca$crit[length(res_rgcca$crit)])

  res_perm <- rgcca_permutation_k(
    res_rgcca2$call, perm = FALSE, par_type = "tau",
    par_value = rep(1, length(blocks))
  )$crit
  crit <- sum(vapply(
    res_rgcca2$crit, function(x) x[length(x)],
    FUN.VALUE = double(1)
  ))
  expect_equal(res_perm, crit)
})

test_that("rgcca_permutation_k gives a smaller criterion than rgcca if perm
          is TRUE", {
  res_perm <- rgcca_permutation_k(
    res_rgcca$call, perm = TRUE, par_type = "tau",
    par_value = rep(1, length(blocks))
  )$crit
  expect_lt(res_perm, res_rgcca$crit[length(res_rgcca$crit)])

  res_perm <- rgcca_permutation_k(
    res_rgcca2$call, perm = TRUE, par_type = "tau",
    par_value = rep(1, length(blocks))
  )$crit
  crit <- sum(vapply(
    res_rgcca2$crit, function(x) x[length(x)],
    FUN.VALUE = double(1)
  ))
  expect_lt(res_perm, crit)
})

test_that("rgcca_permutation_k gives the same criterion as sgcca if perm
          is FALSE", {
  res_perm <- rgcca_permutation_k(
    res_sgcca$call, perm = FALSE, par_type = "sparsity", par_value = c(.75, .8)
  )$crit
  expect_equal(res_perm, res_sgcca$crit[length(res_sgcca$crit)])

  res_perm <- rgcca_permutation_k(
    res_sgcca2$call, perm = FALSE, par_type = "sparsity", par_value = c(.75, .8)
  )$crit
  crit <- sum(vapply(
    res_sgcca2$crit, function(x) x[length(x)],
    FUN.VALUE = double(1)
  ))
  expect_equal(res_perm, crit)
})

test_that("rgcca_permutation_k gives a smaller criterion than sgcca if perm
          is TRUE", {
  res_perm <- rgcca_permutation_k(
    res_sgcca$call, perm = TRUE, par_type = "sparsity", par_value = c(.75, .8)
  )$crit
  expect_lt(res_perm, res_sgcca$crit[length(res_sgcca$crit)])

  res_perm <- rgcca_permutation_k(
    res_sgcca2$call, perm = TRUE, par_type = "sparsity", par_value = c(.75, .8)
  )$crit
  crit <- sum(vapply(
    res_sgcca2$crit, function(x) x[length(x)],
    FUN.VALUE = double(1)
  ))
  expect_lt(res_perm, crit)
})
