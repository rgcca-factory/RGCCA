data("Russett")
blocks <- list(
    agriculture = Russett[, seq(3)],
    industry = Russett[, 4:5],
    politic = Russett[, 6:11] )
rgcca_out <- rgcca(blocks, sparsity = 0.75, type = "sgcca")
boot <- bootstrap(rgcca_out, 2, n_cores = 1)
selected.var <- get_bootstrap(boot, n_cores = 1)

test_that("plot_boot_default", {
    expect_is(plot_bootstrap_2D(boot, n_cores = 1), "ggplot")
})   

rgcca_out <- rgcca(blocks)
boot <- bootstrap(rgcca_out, 2, n_cores = 1)
selected.var <- get_bootstrap(boot, n_cores = 1)

test_that("plot_boot_with_args", {
    expect_is(plot_bootstrap_2D(boot, n_cores = 1), "ggplot")
    expect_is(plot_bootstrap_2D(df_b = selected.var), "ggplot")
    expect_is(plot_bootstrap_1D(boot, n_cores = 1), "ggplot")
    expect_is(plot_bootstrap_2D(df_b = selected.var), "ggplot")
})

rgcca_out <- rgcca(
    blocks,ncomp = c(2,2,2),
    scheme = function(x) x^4, 
    type = "sgcca",
    sparsity = c(.6, .75, .5, 1)
)

boot <- bootstrap(rgcca_out, 2, n_cores = 1)
test_that("plot_boot_object", {
    expect_is(plot(boot, i_block = 4, n_cores = 1), "ggplot")
    expect_is(plot(boot, type = "2D", n_cores = 1), "ggplot")
})
