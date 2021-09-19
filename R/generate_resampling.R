# Product for Matrix (pm) is a generalization of the matricial product %*% for matrices with missing data. Missing data are replaced by 0.
# @param M1  A matrix with n1 lines and p columns
# @param M2  A matrix with p lines and n2 columns
# @param na.rm if TRUE calculates the matricial product only on available data. Else returns NA.
# @return \item{X}{The resulting matrix with n1 lines and n2 columns}
# @title Product for Matrices with missing data (pm)

generate_resampling <- function(rgcca_res, n_boot){
  N                   = NROW(rgcca_res$call$raw[[1]])
  full_idx            = rep(x = seq(N), each = n_boot)
  full_idx            = sample(x = full_idx, size = N*n_boot)
  full_idx            = split(full_idx, ceiling(seq_along(full_idx)/N))
  boot_blocks         = lapply(full_idx, function(idx) lapply(rgcca_res$call$raw,
                                                              function(x){
                                                                y           = x[idx, , drop = FALSE]
                                                                rownames(y) = paste("S",1:length(idx))
                                                                return(y)}))
  boot_column_sd_null = lapply(boot_blocks, function(blocks) remove_null_sd(blocks)$column_sd_null)
  eval_boot_sample    = sapply(boot_column_sd_null, function(x) sum(sapply(x, length)))
  if (sum(eval_boot_sample) == 0){
    summarize_column_sd_null = NULL
  }else{
    summarize_column_sd_null = Reduce("rbind", boot_column_sd_null)
    summarize_column_sd_null = apply(summarize_column_sd_null, 2,
                                     function(x){
                                       vec_x = unlist(x)
                                       dup_x = duplicated(vec_x)
                                       return(vec_x[which(!dup_x)])
                                     })
    boot_blocks              = lapply(boot_blocks,
                                      function(blocks)
                                        remove_null_sd(blocks, summarize_column_sd_null)$list_m)
    warning(paste("Variables : ",
                  paste(names(Reduce("c", summarize_column_sd_null)), collapse = " - "),
                  "appear to be of null variance in some bootstrap sample.\n",
                  "==> Consequently, they were removed from their corresponding blocks.\n",
                  "==> RGCCA is going to be run again without these variables."))
  }
  return(list(boot_blocks = boot_blocks, full_idx = full_idx, summarize_column_sd_null = summarize_column_sd_null))
}

##Future
# id_boot        = sample(NROW(blocks[[1]]), replace = TRUE)
# boot_blocks    = lapply(blocks, function(x) x[id_boot, , drop = FALSE])
# is_null_sd_var = Reduce("|",
#                         unlist(sapply(blocks, function(x)
#                           apply(x, 2, function(x) sd(x, na.rm = TRUE)) == 0)))
# if (is_null_sd_var){
#   browser()
#   probs = sapply(blocks, function(block) apply(block, 2, function(var) {
#     occurences = table(var, useNA = "ifany")/length(var)
#     new_idx    = match(as.character(var), names(occurences))
#     new_var    = as.matrix(occurences[new_idx])
#     new_var    = (1-new_var)/sum(1-new_var)
#     return(new_var)
#   }))
#   probs = apply(Reduce('cbind', probs), 1, max)/sum(apply(Reduce('cbind', probs), 1, max))
#   iter = 0
#   while((is_null_sd_var) && (iter <= 5)){
#     id_boot        = sample(NROW(blocks[[1]]), replace = TRUE, prob = probs)
#     boot_blocks    = lapply(blocks, function(x) x[id_boot, , drop = FALSE])
#     is_null_sd_var = Reduce("|",
#                             unlist(sapply(blocks, function(x)
#                               apply(x, 2, function(x) sd(x, na.rm = TRUE)) == 0)))
#     iter           = iter + 1
#   }
#   if (iter > 5){
#     stop_rgcca("Impossible to find a bootstrap sample without variables",
#                "with null standard deviation. Please take care of
#                these variable, especially : ")
#   }
# }
