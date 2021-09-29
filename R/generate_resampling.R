#' Generate n_boot block-wise bootstrapped samples.
#' @inheritParams bootstrap
#' @param balanced A boolean indicating if a balanced bootstrap procedure is
#' performed or not (default is TRUE).
#' @param keep_all_variables A boolean indicating if all variables have to be
#' kept even when some of them have null variance for at least one bootstrap
#' sample (default is FALSE).
#' @return \boot_blocks}{A list of size n_boot containing all the
#' bootstrapped blocks.}
#' @return \item{full_idx}{A list of size n_boot containing the observations
#' kept for each bootstrap sample.}
#' @return \item{summarize_column_sd_null}{A list of size the number of block
#' containing the variables that were removed from each block for all the
#' bootstrap samples. Variables are removed if they appear to be of null variance
#' in at least one bootstrap sample. If no variable is removed, return NULL.}
#' @title Generate sampled blocks for bootstrap.
#' @keywords internal

generate_resampling <- function(rgcca_res, n_boot, balanced = TRUE, keep_all_variables = FALSE){
  NO_null_sd_var = FALSE
  iter           = 0
  N              = NROW(rgcca_res$call$raw[[1]])
  prob           = rep(1/N, N)
  while (!NO_null_sd_var){
    if (balanced){
      full_idx = rep(x = seq(N), each = n_boot)
      full_idx = sample(x = full_idx, size = N*n_boot, replace = FALSE)
      full_idx = split(full_idx, ceiling(seq_along(full_idx)/N))
    }else{
      full_idx = lapply(seq(n_boot), function(x) sample(x = seq(N), replace = TRUE, prob = prob))
    }

    boot_blocks         = lapply(full_idx, function(idx) lapply(rgcca_res$call$raw,
                                                                function(x){
                                                                  y           = x[idx, , drop = FALSE]
                                                                  rownames(y) = paste("S",1:length(idx))
                                                                  return(y)}))
    boot_column_sd_null = lapply(boot_blocks, function(blocks) remove_null_sd(blocks)$column_sd_null)
    eval_boot_sample    = sapply(boot_column_sd_null, function(x) sum(sapply(x, length)))
    NO_null_sd_var      = (sum(eval_boot_sample) == 0)
    if (NO_null_sd_var){
      summarize_column_sd_null = NULL
    }else{
      if (!keep_all_variables){
        NO_null_sd_var           = TRUE
        summarize_column_sd_null = Reduce("rbind", boot_column_sd_null)
        summarize_column_sd_null = apply(as.data.frame(summarize_column_sd_null), 2,
                                         function(x) unique(x)[1:2][[2]])
        boot_blocks              = lapply(boot_blocks,
                                          function(blocks)
                                            remove_null_sd(blocks, summarize_column_sd_null)$list_m)
        warning(paste("Variables : ",
                      paste(names(Reduce("c", summarize_column_sd_null)), collapse = " - "),
                      "appear to be of null variance in some bootstrap samples",
                      "and thus were removed from all samples. \n",
                      "==> RGCCA is run again without these variables."))
      }else{
        if (iter > 5){
          summarize_column_sd_null = Reduce("rbind", boot_column_sd_null)
          summarize_column_sd_null = apply(as.data.frame(summarize_column_sd_null), 2,
                                           function(x) unique(x)[1:2][[2]])
          error_message = paste("Impossible to define all bootstrap samples",
                                "without variables with null variance. Please",
                                "consider removing these variables: ",
                                paste(names(Reduce("c", summarize_column_sd_null)), collapse = " - "))
          if (balanced){
            error_message = paste0(error_message,
                                   ". Please, consider unbalanced bootstrap by",
                                   " setting 'balanced' to FALSE.")
          }
          stop_rgcca(error_message)
        }
        if (!balanced){
          if (iter == 0){
            prob = sapply(rgcca_res$call$raw,
                           function(block) apply(block, 2, function(var) {
                             occurences = table(var, useNA = "ifany")/length(var)
                             new_idx    = match(as.character(var), names(occurences))
                             new_var    = as.matrix(occurences[new_idx])
                             new_var    = (1-new_var)/sum(1-new_var)
                             return(new_var)
                           }))
            prob = apply(Reduce('cbind', prob), 1, max)/sum(apply(Reduce('cbind', prob), 1, max))
          }
        }
        iter = iter + 1
      }
    }
  }
  return(list(boot_blocks = boot_blocks, full_idx = full_idx, summarize_column_sd_null = summarize_column_sd_null))
}
