 library(RGCCA)
 data("Russett")
 blocks = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
     politic = Russett[, 6:11] )
 rgcca_out = rgcca(blocks)
 bootstrap_k(rgcca_out)
 bootstrap_k(rgcca_out, lapply(blocks, scale), superblock = FALSE)
