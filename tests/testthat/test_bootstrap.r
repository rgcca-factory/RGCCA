#'# test bootstrap

#'''

 data("Russett")
 blocks = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
     politic = Russett[, 6:11] )
 rgcca_out = rgcca(blocks)
 bootstrap(rgcca_out, n_boot = 2, n_cores = 1)
 bootstrap(rgcca_out, n_boot = 2, n_cores = 1, blocks = lapply(blocks, scale),
  superblock = FALSE)
