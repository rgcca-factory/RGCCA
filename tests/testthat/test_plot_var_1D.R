
# weights = lapply(seq(3), function(x) matrix(runif(7*2), 7, 2))
# for (i in seq(3))
#row.names(weights[[i]]) <- paste0(letters[i],
#      letters[seq(NROW(weights[[i]]))])
 # TODO : solve the problem
 # weights[[4]] = Reduce(rbind, weights)
 # rgcca_out = list(a = weights)
 # names(rgcca_out$a) = LETTERS[seq(3)]
 # rgcca_out$blocks = lapply(rgcca_out$a, t)
 # rgcca_out$superblock = TRUE
 # # With the 1rst component of the superblock
 # plot_var_1D(rgcca_out, 1, type = "weight")
 # With the 2nd component of the 1rst block by selecting the ten higher weights
 #plot_var_1D(rgcca_out, 2, 10, 1, type = "weight")
 library(RGCCA)
 data("Russett")
 blocks = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
     politic = Russett[, 6:11] )
 rgcca_out = rgcca.analyze(blocks)
 plot_var_1D(rgcca_out, collapse = TRUE)
