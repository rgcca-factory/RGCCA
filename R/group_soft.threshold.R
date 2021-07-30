# The function soft.threshold() soft-thresholds a vector such that the L1-norm constraint is satisfied.
# @param x A numeric vector.
# @param sumabs A numeric constraint on x's L1 norm.
#
# @return Returns a vector resulting from the soft thresholding of \eqn{x} given sumabs
# @keywords manip
#
group_soft.threshold <- function(x,sumabs=1, group_sparsity, var_order){
  x_sub           = sapply(group_sparsity, function(group) norm(x[group], type = "2"))
  proj            = proj_l1_l2(argu = x_sub, a = sumabs)
  x_sub_projected = soft(x = x_sub, d = proj$lambda)
  x_projected     = lapply(1:length(group_sparsity), function(idx_group)
    as.matrix(x_sub_projected[idx_group]*x[group_sparsity[[idx_group]]]/norm(x[group_sparsity[[idx_group]]], type = "2")))
  x_projected     = Reduce("rbind", x_projected)
  x_projected     = x_projected[var_order]
  if (proj$l2_SAT){
    return(x_projected/norm(x_projected, type = "2"))
  }else{
    return(x_projected)
  }
}
