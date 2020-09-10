data("Russett")
blocks <- list(
    agriculture = Russett[, seq(3)],
    industry = Russett[, 4:5],
    politic = Russett[, 6:11])

blocks3=lapply(blocks,scale)
blocks3=lapply(blocks3,function(x) return(x/sqrt(ncol(x))))
blocks2=RGCCA:::scaling(blocks,scale=TRUE,scale_block=TRUE,bias=FALSE)
test_that("scaling_default_1", {
    expect_true(sum(abs(blocks3[[2]]-blocks2[[2]]))<3e-15)})

