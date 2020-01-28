#'# check_get_ctr test

#'''
data("Russett")
 blocks = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
     politic = Russett[, 6:11] )
 rgcca_out = rgcca(blocks, ncomp = c(3,2,4))
 get_ctr(rgcca_out)
 # On the first block and with weights
 get_ctr(rgcca_out, 2, 1, i_block = 1, type = "weight")
 # With 3 components and on the variables of two blocks
 superblocks <- rep(list(Reduce(cbind, c(blocks[1], blocks[3]))), 2)
 names(superblocks) <- names(blocks)[c(1, 3)]
 rgcca_out = rgcca(blocks[c(1,3)], ncomp = c(3,4))

# get_ctr(rgcca_out, compz = 3, i_block = 1, type = "cor", collapse = TRUE)
# rgcca_out$blocks = superblocks
# get_ctr(rgcca_out, compx=2, compy=1, compz=3, i_block=1, "weight", collapse=TRUE)