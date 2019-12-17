#  data("Russett")
#  blocks = list(agriculture = Russett[, seq(3)],
#      politic = Russett[, 6:11] )
#  rgcca_out = rgcca.analyze(blocks, ncomp = rep(3, 2))
# df = get_comp(rgcca_out, compz = 3)
# plot3D(df, rgcca_out, i_block = 2)
#  plot3D(df, rgcca_out, i_block = 2, text = FALSE)
#  response = factor( apply(Russett[, 9:11], 1, which.max),
#                    labels = colnames(Russett)[9:11] )
# response = blocks[[2]][, 1]
#  names(response) = row.names(blocks[[2]])
#  df = get_comp(rgcca_out, response, compz = 3)
#  plot3D(df, rgcca_out, i_block = 2, text = FALSE)
#  plot3D(df, rgcca_out, i_block = 2)
#  df = get_ctr2(rgcca_out, compz = 3, i_block = 1, collapse = TRUE)
#  plot3D(df, rgcca_out, i_block = 2, type = "var")
