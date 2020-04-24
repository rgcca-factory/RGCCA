set.seed(1)

data("Russett")
blocks <- list(
    agriculture = Russett[, seq(3)],
    industry = Russett[, 4:5],
    politic = Russett[, 6:11])

# Superblock=TRUE
#-------------------------
test_that("set_rgcca_equal_to_rgcca", {
    rgcca_out <- rgcca(blocks, response = 2,superblock=TRUE)
    res <- set_rgcca(
        rgcca_out,
        inds =  .Machine$integer.max,
        blocks = blocks,
        response = 2,
        tol = 1E-8,
        superblock=TRUE
    )
    expect_identical(round(res$Y[[1]], 7), round(rgcca_out$Y[[1]], 7))
}) # TODO ! !!

rgcca_out <- rgcca(blocks, response = 2)
res1 <- set_rgcca(
    rgcca_out,
    inds =  .Machine$integer.max,
    blocks = blocks,
    response = 2,
    tol = 1E-8
)
res2<- set_rgcca(
    rgcca_out,
    inds =  .Machine$integer.max,
    tol = 1E-8
)
cor(res1$Y[[1]][,1],res2$Y[[1]][,1])
res1$Y[[2]][,2]==res2$Y[[2]][,2]

# Superblock=FALSE
#------------------
# Without scaling
rgcca_out <- rgcca(blocks, response = 2,superblock=FALSE,scale=FALSE,sameBlockWeight = FALSE)
res <- set_rgcca(
    rgcca_out,
    inds =  .Machine$integer.max,
    blocks = blocks,
    response = 2,
    tol = 1E-8,
    superblock=FALSE,
    ,scale=FALSE,sameBlockWeight = FALSE
)
# scale = TRUE
round(rgcca_out$Y[[1]][,1],digits=7)==round(res$Y[[1]][,1],digits=7)
rgcca_out <- rgcca(blocks, response = 2,superblock=FALSE,scale=TRUE,sameBlockWeight = FALSE)
res <- set_rgcca(
    rgcca_out,
    inds =  .Machine$integer.max,
    blocks = blocks,
    response = 2,
    tol = 1E-8,
    superblock=FALSE,
    ,scale=TRUE,sameBlockWeight = FALSE
)
round(rgcca_out$Y[[1]][,1],digits=7)==round(res$Y[[1]][,1],digits=7)
