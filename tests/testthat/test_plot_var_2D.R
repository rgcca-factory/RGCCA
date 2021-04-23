
# setMatrix = function(nrow, ncol, iter = 3) lapply(seq(iter),
#    function(x) matrix(runif(nrow * ncol), nrow, ncol))
# blocks = setMatrix(10, 5)
#  blocks[[4]] = Reduce(cbind, blocks)
# for (i in seq(4)) {
#     colnames(blocks[[i]]) = paste0( LETTERS[i],
#     as.character(seq(NCOL(blocks[[i]]))))
# }
#  coord = setMatrix(10, 2, 4)
#  a = setMatrix(5, 2)
#  a[[4]] = matrix(runif(15 * 2), 15, 2)
#  AVE_X = lapply(seq(4), function(x) runif(2))
#  rgcca_out = list(Y = coord, a = a, AVE = list(AVE_X = AVE_X), blocks = blocks)
#  names(rgcca_out$a) <- LETTERS[seq(4)] -> names(rgcca_out$blocks)
#  # Using a superblock
#  rgcca_out$superblock = TRUE
#  class(rgcca_out) = "rgcca"
#  plot_var_2D(rgcca_out, 1, 2)
#  # Using the first block
#  plot_var_2D(rgcca_out, 1, 2, 1)
 data("Russett")
 blocks = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
    politic = Russett[, 6:11] )
rgcca_out = rgcca(blocks, ncomp = 2)
 plot_var_2D(rgcca_out, collapse = TRUE)
 plot_var_2D(rgcca_out)

 RussettWithNA=Russett
 RussettWithNA[1:2,1:3]=NA
 RussettWithNA[3,4:5]=NA
 RussettWithNA[3,1]=NA
 blocksNA = list(agriculture = RussettWithNA[, seq(3)], industry = RussettWithNA[, 4:5],
                 politic = RussettWithNA[, 6:11] )
 resRGCCANA1=rgcca(blocksNA,NA_method="complete", ncomp = 2)
 plot_var_2D(resRGCCANA1, collapse = TRUE)
