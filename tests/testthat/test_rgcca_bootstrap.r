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
  expect_is(res, "bootstrap")
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
select_var <- dplyr::filter(
  resBootstrap$stats, var == "demostab", type == "weights", comp == 1
)
test_that("test_rgcca_bootstrap_na_values", {
  expect_equal(
    select_var$mean,
    mean(dplyr::filter(
      resBootstrap$bootstrap, var == "demostab", type == "weights", comp == 1
    )$value)
  )
  expect_equal(
    select_var$estimate, resRGCCA$a[[3]]["demostab", 1]
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

##############################################
# Test on the risk of having null variance   #
# variables in at least one bootstrap sample #
##############################################
# Here, the variable `rent` is trapped and should be detected in
# `generate_resampling` which is going to raise a warning. Then it
# should be removed by `bootstrap`.
blocks <- list(
  agriculture = Russett[, seq(3)],
  industry    = Russett[, 4:5],
  politic     = Russett[, 6:11]
)

ncomp <- 1
# Rent is trapped.
blocks$agriculture$rent <- 0
blocks$agriculture$rent[1:4] <- 1
rgcca_out <- rgcca(blocks, ncomp = ncomp)

set.seed(8882)
test_that(
  "rgcca_bootstrap_removed_variable_1",
  expect_warning(
    rgcca_bootstrap(rgcca_out, n_boot = 4, n_cores = 1, balanced = TRUE),
    paste0(
      "Variables:  rent appear to be of null ",
      "variance in some bootstrap samples and thus ",
      "were removed from all samples. \n",
      " ==> RGCCA is run again without these variables."
    )
  )
)
# Exact same situation where we check in the output that `rent` was removed.
set.seed(8882)

test_that("rgcca_bootstrap_removed_variable_2", {
  expect_warning(
    boot_out <- rgcca_bootstrap(
      rgcca_out, n_boot = 4, n_cores = 1, balanced = TRUE
    ),
    paste0(
      "Variables:  rent appear to be of null variance in some bootstrap ",
      "samples and thus were removed from all samples."
    ),
    fixed = TRUE
  )
  expect_false("rent" %in% boot_out$bootstrap$var)
  expect_false("rent" %in% colnames(boot_out$rgcca$blocks$agriculture))
  expect_false("rent" %in% colnames(boot_out$rgcca$call$blocks$agriculture))
})

# Same situation, but this time, it is specifically ask that all variables are
# kept. It is thus checked that `rent` is still there.
set.seed(8882)
boot_out <- rgcca_bootstrap(rgcca_out,
  n_boot = 4, n_cores = 1,
  keep_all_variables = TRUE, balanced = TRUE
)

test_that("rgcca_bootstrap_keep_all_variables", {
  expect_true("rent" %in% boot_out$bootstrap$var)
  expect_true("rent" %in% colnames(boot_out$rgcca$blocks$agriculture))
  expect_true("rent" %in% colnames(boot_out$rgcca$call$blocks$agriculture))
})
