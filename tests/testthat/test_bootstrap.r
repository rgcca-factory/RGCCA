data("Russett")
blocks <- list(
  agriculture = Russett[, seq(3)],
  industry = Russett[, 4:5],
  politic = Russett[, 6:11]
)
p <- vapply(blocks, NCOL, FUN.VALUE = 1L)
ncomp <- 1
rgcca_out <- rgcca(blocks, ncomp = 1)
boot <- bootstrap(rgcca_out, n_boot = 4, n_cores = 1)

test_that("bootstrap_default_1", {
  expect_equal(length(boot), 2)
  expect_equal(length(boot$bootstrap), 2)
  boot1 <- boot$bootstrap[[1]][[1]]
  expect_is(boot, "bootstrap")
  expect_is(boot$rgcca, "rgcca")
  expect_is(boot1, "list")
  expect_is(boot1[[1]], "matrix")
  expect_true(all(vapply(boot1, NCOL, FUN.VALUE = 1L) == 4))
  expect_identical(vapply(boot1, NROW, FUN.VALUE = 1L), p)
})

### Case ncomp = 2
rgcca_out <- rgcca(blocks, ncomp = 2)
boot <- bootstrap(rgcca_out, n_boot = 2, n_cores = 1)

test_that("bootstrap_default", {
  expect_equal(length(boot), 2)
  expect_equal(length(boot$bootstrap), 2)
  boot1 <- boot$bootstrap[[1]][[1]]
  expect_is(boot, "bootstrap")
  expect_is(boot$rgcca, "rgcca")
  expect_is(boot1, "list")
  expect_is(boot1[[1]], "matrix")
  expect_true(all(vapply(boot1, NCOL, FUN.VALUE = 1L) == 2))
  expect_identical(vapply(boot1, NROW, FUN.VALUE = 1L), p)
})

### Case sgcca
rgcca_out <- rgcca(blocks,
  ncomp = 2, method = "sgcca",
  sparsity = c(0.8, 0.9, 0.7), superblock = TRUE
)
boot <- bootstrap(rgcca_out, n_boot = 2, n_cores = 1)

test_that("bootstrap_default", {
  expect_equal(length(boot), 2)
  expect_equal(length(boot$bootstrap), 2)
  boot1 <- boot$bootstrap[[1]][[1]]
  expect_is(boot, "bootstrap")
  expect_is(boot$rgcca, "rgcca")
  expect_is(boot1, "list")
  expect_is(boot1[[1]], "matrix")
  expect_true(all(vapply(boot1, NCOL, FUN.VALUE = 1L) == 2))
})

blocks[[1]][1:3, 1] <- NA
blocks[[1]][4, ] <- NA
resRGCCA <- rgcca(blocks, ncomp = c(2, 2, 2), superblock = FALSE)
set.seed(seed = 18)
resBootstrap <- bootstrap(rgcca = resRGCCA, n_boot = 2, n_cores = 1)
select_var <- get_bootstrap(resBootstrap, display_order = TRUE)
test_that("test_bootstrap_na_values", {
  expect_equal(
    select_var["demostab", "mean"],
    mean(c(resBootstrap$bootstrap[[1]][[1]][["politic"]]["demostab", ]))
  )
  expect_true(
    select_var["demostab", "estimate"] == resRGCCA$a[[3]]["demostab", 1]
  )
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
  "bootstrap_removed_variable_1",
  expect_warning(
    bootstrap(rgcca_out, n_boot = 4, n_cores = 1, balanced = T),
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

test_that("bootstrap_removed_variable_2", {
  expect_warning(
    boot_out <- bootstrap(rgcca_out, n_boot = 4, n_cores = 1, balanced = T),
    paste0(
      "Variables:  rent appear to be of null variance in some bootstrap ",
      "samples and thus were removed from all samples."
    ),
    fixed = TRUE
  )
  expect_false("rent" %in% rownames(boot_out$bootstrap$W[[1]]$agriculture))
  expect_false("rent" %in% rownames(boot_out$bootstrap$L[[1]]$agriculture))
  expect_false("rent" %in% colnames(boot_out$rgcca$call$blocks$agriculture))
  expect_false("rent" %in% colnames(boot_out$rgcca$call$raw$agriculture))
})

# Same situation, but this time, it is specifically ask that all variables are
# kept. It is thus checked that `rent` is still there.
set.seed(8882)
boot_out <- bootstrap(rgcca_out,
  n_boot = 4, n_cores = 1,
  keep_all_variables = T, balanced = T
)

test_that("bootstrap_keep_all_variables", {
  expect_true("rent" %in% rownames(boot_out$bootstrap$W[[1]]$agriculture))
  expect_true("rent" %in% rownames(boot_out$bootstrap$L[[1]]$agriculture))
  expect_true("rent" %in% colnames(boot_out$rgcca$call$blocks$agriculture))
  expect_true("rent" %in% colnames(boot_out$rgcca$call$raw$agriculture))
})
