data("Russett")
blocks <- list(
  agriculture = Russett[, seq(3)],
  industry = Russett[, 4:5],
  politic = Russett[, 6:11]
)
p <- vapply(blocks, NCOL, FUN.VALUE = 1L)

test_structure <- function(res, n_boot, ncomp, p) {
  expect_equal(length(res), 4)
  expect_equal(dim(res$bootstrap), c(2 * n_boot * ncomp * sum(p), 6))
  expect_is(res, "rgcca_bootstrap")
  expect_is(res$rgcca, "rgcca")
}

test_that("rgcca_bootstrap_default_1", {
  rgcca_out <- rgcca(blocks, ncomp = 1)
  boot <- rgcca_bootstrap(rgcca_out, n_boot = 4, n_cores = 1)
  test_structure(boot, 4, 1, p)
})

### Case ncomp = 2
test_that("rgcca_bootstrap_default", {
  rgcca_out <- rgcca(blocks, ncomp = 2, tau = "optimal")
  boot <- rgcca_bootstrap(rgcca_out, n_boot = 2, n_cores = 1)
  test_structure(boot, 2, 2, p)
})

### Case sgcca
test_that("rgcca_bootstrap_default", {
  rgcca_out <- rgcca(blocks,
                     ncomp = 2, method = "sgcca",
                     sparsity = c(0.8, 0.9, 0.7), superblock = TRUE
  )
  q <- vapply(rgcca_out$a[-4], function(x) {
    length(unique(which(x != 0, arr.ind = TRUE)[, 1]))
  }, FUN.VALUE = 1L)
  q <- c(q, sum(q))
  boot <- rgcca_bootstrap(rgcca_out, n_boot = 2, n_cores = 1)
  test_structure(boot, 2, 2, q)
})

blocks[[1]][1:3, 1] <- NA
blocks[[1]][4, ] <- NA
resRGCCA <- rgcca(blocks, ncomp = c(2, 2, 2), superblock = FALSE)
set.seed(seed = 18)
resBootstrap <- rgcca_bootstrap(rgcca = resRGCCA, n_boot = 2, n_cores = 1)
select_var <- subset(
  resBootstrap$stats, var == "demostab" & type == "weights" & comp == 1
)
test_that("test_rgcca_bootstrap_na_values", {
  expect_equal(
    select_var$mean,
    mean(subset(
      resBootstrap$bootstrap, var == "demostab" & type == "weights" & comp == 1
    )$value)
  )
  expect_equal(
    select_var$estimate, unname(resRGCCA$a[[3]]["demostab", 1])
  )
})

### Case qualitative block
test_that("rgcca_bootstrap_classif", {
  blocks_classif <- list(
    agriculture = Russett[, seq(3)],
    industry = Russett[, 4:5],
    politic = as.factor(apply(Russett[, 9:11], 1, which.max))
  )
  rgcca_out <- rgcca(blocks_classif, response = 3)
  p <- vapply(rgcca_out$blocks, NCOL, FUN.VALUE = 1L)
  boot <- rgcca_bootstrap(rgcca_out, n_boot = 4, n_cores = 1)
  test_structure(boot, 4, 1, p)
})
