#'# get_block_var
#'''
#' library(RGCCA)
 data("Russett")
 blocks = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
     politic = Russett[, 6:11] )
 rgcca_out = rgcca(blocks)
 response = factor( apply(Russett[, 9:11], 1, which.max),
                   labels = colnames(Russett)[9:11] )
 get_comp(rgcca_out, as.matrix(response))
 response = as.matrix(runif(NROW(blocks[[1]])))
 row.names(response) = row.names(blocks[[1]])
 get_comp(rgcca_out, response)
