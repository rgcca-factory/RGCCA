#'# test bootstrap

#'''

data("Russett")
block <- list(
    agriculture = Russett[, seq(3)],
    industry = Russett[, 4:5],
    politic = Russett[, 6:11] )
rgcca_out <- rgcca(blocks)
boot <- bootstrap(rgcca_out, n_boot = 2, n_cores = 1)

test_that("bootstrap_default", {
    expect_equal(length(boot), 2)
    expect_equal(length(boot$bootstrap), 2)
    boot1 <- boot$bootstrap[[1]]
    expect_is(boot, "bootstrap")
    expect_is(boot$rgcca, "rgcca")
    expect_is(boot1, "list")
    expect_is(boot1[[1]], "matrix")
    expect_true(all(sapply(boot1, NCOL) == 2))
    expect_identical(sapply(boot1, NROW), sapply(rgcca_out$call$blocks, NCOL))
})

test_that("bootstrap_default", {
    select_var <- get_bootstrap(boot, n_cores = 1)
    expect_is(select_var, "df_bootstrap")
    expect_is(select_var, "data.frame")
    expect_identical(NROW(select_var), NCOL(rgcca_out$call$blocks[[4]]))
})

test_that("bootstrap_with_args", {
    expect_is(
        bootstrap(
            rgcca_out, 
            n_boot = 2, 
            n_cores = 1, 
            blocks = lapply(blocks, scale),
            superblock = FALSE),
        "bootstrap")
})

blocks[[1]][1:3, 1] <- NA
blocks[[1]][4,] <- NA
resRGCCA <- rgcca(blocks, ncomp = c(2,2,2))
set.seed(seed = 18)
resBootstrap <- bootstrap( rgcca=resRGCCA, n_boot = 2, n_cores = 1)
select_var <- get_bootstrap(resBootstrap, n_cores = 1)
plot_bootstrap_1D(df_b = select_var)

test_that("test_bootstrap_na_values", {
    expect_equal(
        select_var["labo", 1],
        mean(c(resBootstrap$bootstrap[[1]][[4]]["labo", 1], resBootstrap$bootstrap[[2]][[4]]["labo", 1]))
    )
    expect_true(select_var["labo", 2] == resRGCCA$a[[4]]["labo", 1])
})
