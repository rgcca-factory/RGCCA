# The function soft.threshold() soft-thresholds a vector such that the L1-norm constraint is satisfied.
# @param x A numeric vector.
# @param sumabs A numeric constraint on x's L1 norm.
#
# @return Returns a vector resulting from the soft thresholding of \eqn{x} given sumabs
# @keywords manip
#
group_soft.threshold <- function(x,sumabs=1, group_sparsity){
  x_sub           = sapply(group_sparsity, function(group) norm(x[group], type = "2"))
  x_sub_projected = soft(x_sub, proj_l1_l2(x_sub,sumabs)$lambda)
  x_projected     = lapply(1:length(group_sparsity), function(idx_group) 
    as.matrix(x_sub_projected[idx_group]*x[group_sparsity[[idx_group]]]/norm(x[group_sparsity[[idx_group]]], type = "2")))
  x_projected     = Reduce("rbind", x_projected)
  return(x_projected)
}
