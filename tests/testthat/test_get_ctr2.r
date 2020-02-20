#'# check_get_ctr test

#'''
#' library(RGCCA)
 data("Russett")
 blocks = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
     politic = Russett[, 6:11] )
 rgcca.res = rgcca(blocks, ncomp = c(3, 2, 4))
 get_ctr2(rgcca.res, i_block = 2)
 blocks = blocks[c(1,3)]
 rgcca.res = rgcca(blocks, ncomp = c(3,4))
 get_ctr2(rgcca.res, compz = 3, i_block = 1, collapse = TRUE)
 get_ctr2(rgcca.res, 1, 2, 3, 1, "weight", collapse = TRUE, n_mark = 5)
 get_ctr2(rgcca.res, collapse = TRUE)