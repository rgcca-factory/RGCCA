set.seed(1)

data("Russett")
blocks <- list(
    agriculture = Russett[, seq(3)],
    industry = Russett[, 4:5],
    politic = Russett[, 6:11])

# Superblock=TRUE
#-------------------------
rgcca_out <- rgcca(blocks, response = 2,superblock=TRUE)
res <- set_rgcca(
    rgcca_out,
    inds =  .Machine$integer.max,
    blocks = blocks,
    response = 2,
    tol = 1E-8,
    superblock=TRUE
)
round(rgcca_out$Y[[1]][,1],digits=7)==round(res$Y[[1]][,1],digits=7)

test_that("set_rgcca_equal_to_rgcca_sb_t", {
    expect_identical(round(res$Y[[1]], 7), round(rgcca_out$Y[[1]], 7))
}) # TODO ! !!


# Superblock=FALSE: everything is ok
#-------------------------------------
rgcca_out <- rgcca(blocks, response = 2,superblock=FALSE,scale=FALSE,sameBlockWeight = FALSE)
res <- set_rgcca(
    rgcca_out,
    inds =  .Machine$integer.max,
    blocks = blocks,
    response = 2,
    tol = 1E-8,
    superblock=FALSE    ,scale=FALSE,sameBlockWeight = FALSE
)
round(rgcca_out$Y[[1]][,1],digits=7)==round(res$Y[[1]][,1],digits=7)

test_that("set_rgcca_equal_to_rgcca_s_f_sbw_f_sb_f", {
    expect_identical(round(rgcca_out$Y[[1]][,1],digits=7),round(res$Y[[1]][,1],digits=7))
}) 

rgcca_out <- rgcca(blocks, response = 2,superblock=FALSE,scale=TRUE,sameBlockWeight = FALSE)
res <- set_rgcca(
    rgcca_out,
    inds =  .Machine$integer.max,
    blocks = blocks,
    response = 2,
    tol = 1E-8
)
round(rgcca_out$Y[[1]][,1],digits=7)==round(res$Y[[1]][,1],digits=7)
test_that("set_rgcca_equal_to_rgcca_s_t_sbw_f_sb_f", {
    expect_identical(round(rgcca_out$Y[[1]][,1],digits=7),round(res$Y[[1]][,1],digits=7))
}) 
# scale = TRUE & sameBlockWeight=TRUE & superblock=FALSE
rgcca_out <- rgcca(blocks, response = 2,superblock=FALSE,scale=TRUE,sameBlockWeight = TRUE)
res <- set_rgcca(
    rgcca_out,
    inds =  .Machine$integer.max,
    blocks = blocks,
    response = 2,
    tol = 1E-8
    #,
    #superblock=FALSE
    #scale=TRUE,
    #sameBlockWeight = TRUE
)
round(rgcca_out$Y[[1]][,1],digits=7)==round(res$Y[[1]][,1],digits=7)
test_that("set_rgcca_equal_to_rgcca_s_t_sbw_t_sb_f", {
    expect_identical(round(rgcca_out$Y[[1]][,1],digits=7),round(res$Y[[1]][,1],digits=7))
}) 

res <- set_rgcca(
    rgcca_out,
    inds =  .Machine$integer.max,
    blocks = blocks,
    response = 2,
    tol = 1E-8,
    tau=c(0.8,0.9,1)
)
names(blocks)
res$call$tau


