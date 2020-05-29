get_sign_cor = function(i_block, collapse = F, type = "names") { 
  p.vals <- order_df(get_ctr(rgcca_res2, i_block = i_block, type = "cor.test", collapse = collapse), 1)
  i_sign <- which(p.vals[, 1] < .05 / NROW(p.vals))
  if (type == "names")
    rownames(p.vals)[i_sign]
  else
    p.vals[i_sign, ]
}