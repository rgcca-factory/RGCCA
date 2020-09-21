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
rgcca_out <- rgcca(blocks, response = 2,superblock=FALSE,scale=FALSE,scale_block = FALSE)
res <- set_rgcca(
    rgcca_out,
    inds =  .Machine$integer.max,
    blocks = blocks,
    response = 2,
    tol = 1E-8,
    superblock=FALSE    ,scale=FALSE,scale_block = FALSE
)
round(rgcca_out$Y[[1]][,1],digits=7)==round(res$Y[[1]][,1],digits=7)

test_that("set_rgcca_equal_to_rgcca_s_f_sbw_f_sb_f", {
    expect_identical(round(rgcca_out$Y[[1]][,1],digits=7),round(res$Y[[1]][,1],digits=7))
}) 

rgcca_out <- rgcca(blocks, response = 2,superblock=FALSE,scale=TRUE,scale_block = FALSE)
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
# scale = TRUE & scale_block=TRUE & superblock=FALSE
rgcca_out <- rgcca(blocks, response = 2,superblock=FALSE,scale=TRUE,scale_block = TRUE)
res <- set_rgcca(
    rgcca_out,
    inds =  .Machine$integer.max,
    blocks = blocks,
    response = 2,
    tol = 1E-8
    #,
    #superblock=FALSE
    #scale=TRUE,
    #scale_block = TRUE
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


# Remonving one  

#checking set_rgcca without default
rgcca_out <- rgcca(blocks, response = 1,superblock=FALSE,ncomp=1,scale=TRUE,scale_block=TRUE,tol=1e-8)
rgcca_set_1 <- set_rgcca(rgcca_out,tol=1e-8)
test_that("set_rgcca_identical_for_ind0", {
    expect_identical(all.equal(rgcca_out,rgcca_set_1),TRUE)
}) 
    

blocks_2=lapply(blocks,function(x){return(x[-1,])});
blocks_1=lapply(blocks,function(x){return(x[1,])});
names(blocks_2)=names(blocks_1)=names(blocks)
rgcca_out_2 <- rgcca(blocks_2, response = 1,superblock=FALSE,ncomp=1,scale=TRUE,scale_block=TRUE,tol=1e-8)
rgcca_set_2 <- set_rgcca(rgcca_out,inds=1,tol=1e-8)

test_that("set_rgcca_identical_for_ind1", {
    expect_identical(all.equal(rgcca_out_2,rgcca_set_2),TRUE)
}) 


# set_rgcca for bootstrap
data("Russett")
blocks <- list(
    agriculture = Russett[, seq(3)],
    industry = Russett[, 4:5],
    politic = Russett[, 6:11])
rgcca_out <- rgcca(blocks, response = 1,superblock=FALSE,ncomp=1,scale=TRUE,scale_block=TRUE,tol=1e-8)
rgcca_set_boot <- set_rgcca(rgcca_out,inds=1,tol=1e-8,boot=TRUE)

test_that("for_bootstrap_same_nrows", {
    expect_identical(
    lapply(rgcca_set_boot$call$blocks,dim)[[1]][1]==47
    ,TRUE)})
