data("Russett")
blocks <- list(
    agriculture = Russett[, seq(3)],
    industry = Russett[, 4:5],
    politic = Russett[, 6:11] )
ncomp=1
rgcca_out <- rgcca(blocks,ncomp=1)
boot <- bootstrap(rgcca_out, n_boot = 4, n_cores = 1)

#init=Sys.time();boot <- bootstrap(rgcca_out, n_boot = 1000, n_cores = 1);Sys.time()-init
#init=Sys.time();boot <- bootstrap(rgcca_out, n_boot = 1000, para=FALSE);Sys.time()-init
#init=Sys.time();boot <- bootstrap(rgcca_out, n_boot = 1000, para=TRUE);Sys.time()-init


test_that("bootstrap_default_1", {
    expect_equal(length(boot), 2)
    expect_equal(length(boot$bootstrap), 2)
    boot1 <- boot$bootstrap[[1]][[1]]
    expect_is(boot, "bootstrap")
    expect_is(boot$rgcca, "rgcca")
    expect_is(boot1, "list")
    expect_is(boot1[[1]], "matrix")
    expect_true(all(sapply(boot1, NCOL) == 4))
    expect_identical(sapply(boot1, NROW), sapply(rgcca_out$call$blocks, NCOL))
})

# ncomp=2
rgcca_out <- rgcca(blocks,ncomp=2)
boot <- bootstrap(rgcca_out, n_boot = 2, n_cores = 1)

test_that("bootstrap_default", {
    expect_equal(length(boot), 2)
    expect_equal(length(boot$bootstrap), 2)
    boot1 <- boot$bootstrap[[1]][[1]]
    expect_is(boot, "bootstrap")
    expect_is(boot$rgcca, "rgcca")
    expect_is(boot1, "list")
    expect_is(boot1[[1]], "matrix")
    expect_true(all(sapply(boot1, NCOL) == 2))
    expect_identical(sapply(boot1, NROW), sapply(rgcca_out$call$blocks, NCOL))
})

test_that("bootstrap_default", {
    select_var <- get_bootstrap(boot)
    expect_is(select_var, "df_bootstrap")
    expect_is(select_var, "data.frame")
    expect_identical(NROW(select_var), NCOL(rgcca_out$call$blocks[[length(rgcca_out$call$blocks)]]))
})

blocks[[1]][1:3, 1] <- NA
blocks[[1]][4,] <- NA
resRGCCA <- rgcca(blocks, ncomp = c(2,2,2),superblock=FALSE)
set.seed(seed = 18)
resBootstrap <- bootstrap( rgcca=resRGCCA, n_boot = 2, n_cores = 1)
select_var <- get_bootstrap(resBootstrap,display_order=TRUE)
plot_bootstrap_1D(df_b = select_var)
#plot(resBootstrap)
test_that("test_bootstrap_na_values", {
    expect_equal(
        select_var["demostab", "mean"],
        mean(c(resBootstrap$bootstrap[[1]][[1]][["politic"]]["demostab", ]))
    )
    expect_true(select_var["demostab", "estimate"] == resRGCCA$a[[3]]["demostab", 1])
})

# bootstrap superblock
data("Russett")
blocks <- list(
    agriculture = Russett[, seq(3)],
    industry    = Russett[, 4:5],
    politic     = Russett[, 6:11] )

ncomp     = 1
rgcca_out = rgcca(blocks,ncomp=2,superblock=TRUE)
boot      = bootstrap(rgcca_out, n_boot = 100, n_cores = 1)
plot(boot)
print(boot)
plot(rgcca_out, type = "corCircle")

#############################################
# Test on the risk of having null variance  #
# variables in at least on bootstrap sample #
#############################################
# By setting `pval` to `1` in `bootstrap`, `generate_resampling` will not detect
# that the trapped `rent` variable is risky. Thus `bootstrap` should find that
# `rent` in one bootstrap sample is of null variance.
blocks <- list(
    agriculture = Russett[, seq(3)],
    industry    = Russett[, 4:5],
    politic     = Russett[, 6:11])

ncomp                        = 1
# Rent is trapped.
blocks$agriculture$rent      = 0
blocks$agriculture$rent[1:4] = 1
rgcca_out                    = rgcca(blocks, ncomp = ncomp)

# Thus, a warning message should be sent by `bootstrap` for each bootstrap sample
# it sees where `rent` is of null variance.
set.seed(8882)
test_that("bootstrap_missing_val_identification_1",
          expect_warning(bootstrap(rgcca_out, n_boot = 4, n_cores = 1, pval = 1, balanced = T),
                         paste0("This bootstrap sample was discarded as ",
                                "variables: rent, were removed from it because ",
                                "of their null variance in this sample.")))
# Exact same situation where we check that at least for a whole bootstrap
# sample (the troublesome one), the weight/loadings of all the variables from
# all the blocks are set to NAs.
set.seed(8882)
boot_out              = bootstrap(rgcca_out, n_boot = 4, n_cores = 1, pval = 1, balanced = T)
isNA_for_whole_sample = lapply(boot_out$bootstrap,
                               function(x)
                                   which(apply(is.na(Reduce("rbind", x[[1]])), 2, all)))
test_that("bootstrap_missing_val_identification_2", {
    expect_true(length(isNA_for_whole_sample$W) != 0)
    expect_true(length(isNA_for_whole_sample$L) != 0)
})

# Same situation, but this time, if `pval` is set to its default value, the
# variable `rent` is detected and removed. Thus a warning is produced.
set.seed(8882)
test_that("bootstrap_removed_variable_1",
          expect_warning(bootstrap(rgcca_out, n_boot = 4, n_cores = 1, balanced = T),
                         paste0("Variables:  rent appear to be of null ",
                                "variance in some bootstrap samples and thus ",
                                "were removed from all samples. \n",
                                " ==> RGCCA is run again without these variables.")))
# Exact same situation where we check in the output that `rent` was removed.
set.seed(8882)
boot_out = bootstrap(rgcca_out, n_boot = 4, n_cores = 1, balanced = T)

test_that("bootstrap_removed_variable_2",{
    expect_false("rent" %in% rownames(boot_out$bootstrap$W[[1]]$agriculture))
    expect_false("rent" %in% rownames(boot_out$bootstrap$L[[1]]$agriculture))
    expect_false("rent" %in% colnames(boot_out$rgcca$call$blocks$agriculture))
    expect_false("rent" %in% colnames(boot_out$rgcca$call$raw$agriculture))
})

# Same situation, but this time, `pval` is set to its default value and it is
# specifically ask that all variables are kept. It is thus checked that `rent`
# is still there.
set.seed(8882)
boot_out = bootstrap(rgcca_out, n_boot = 4, n_cores = 1, keep_all_variables = T, balanced = T)

test_that("bootstrap_keep_all_variables",{
    expect_true("rent" %in% rownames(boot_out$bootstrap$W[[1]]$agriculture))
    expect_true("rent" %in% rownames(boot_out$bootstrap$L[[1]]$agriculture))
    expect_true("rent" %in% colnames(boot_out$rgcca$call$blocks$agriculture))
    expect_true("rent" %in% colnames(boot_out$rgcca$call$raw$agriculture))
})
