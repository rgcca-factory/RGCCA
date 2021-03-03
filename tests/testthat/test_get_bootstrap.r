# Bootstrap with random data
nv1=10
nv2=100
a=matrix(rnorm(100*nv1),100,nv1);rownames(a)=paste("a",1:100);colnames(a)=paste("s",1:nv1)
b=matrix(rnorm(100*nv2),100,nv2);rownames(b)=paste("a",1:100);colnames(b)=paste("b",1:nv2)
blocks=list(a=a,b=b)
rgcca_out=rgcca(blocks)
n_boot=100
boot <- bootstrap(rgcca_out,n_boot=n_boot,n_cores=1)
res=get_bootstrap(boot)
plot(boot,block=1,bars="sd", n_mark=10)
plot(boot,block=1,bars="stderr", n_mark=10)
plot(boot,block=1,bars="quantile", n_mark=10)

# Bootstrap on Russett
data("Russett")
blocks <- list(
    agriculture = Russett[, seq(3)],
    industry = Russett[, 4:5],
    politic = Russett[, 6:11] )
rgcca_out <- rgcca(blocks,ncomp=2)
boot <- bootstrap(rgcca_out, n_boot = 2, n_cores = 1)


test_that("get_bootstrap_default", {
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
    select_var <- get_bootstrap(boot)
    expect_is(select_var, "df_bootstrap")
    expect_is(select_var, "data.frame")
    expect_identical(NROW(select_var), NCOL(rgcca_out$call$blocks[[length(rgcca_out$call$blocks)]]))
})
#
# test_that("bootstrap_with_args", {
#     rgcca_out <- rgcca(blocks, superblock = FALSE,ncomp=2)
#     expect_is(
#         bootstrap(
#             rgcca_out,
#             n_boot = 2,
#             n_cores = 1,
#             blocks = lapply(blocks, scale),
#             superblock = FALSE),
#         "bootstrap")
# })

blocks[[1]][1:3, 1] <- NA
blocks[[1]][4,] <- NA
resRGCCA <- rgcca(blocks,superblock=FALSE,ncomp=2)
set.seed(seed = 18)
resBootstrap <- bootstrap( rgcca=resRGCCA, n_boot = 2, n_cores = 1)
select_var <- get_bootstrap(resBootstrap ,display_order=TRUE)
plot_bootstrap_1D(df_b = select_var)

test_that("test_bootstrap_na_values", {
    expect_equal(
        select_var["demostab", 1],
        mean(c(resBootstrap$bootstrap[[1]][["politic"]]["demostab", ]))
    )
    expect_true(select_var["demostab", "estimate"] == resRGCCA$a[[3]]["demostab", 1])
})
