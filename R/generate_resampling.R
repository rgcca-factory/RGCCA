#' Generate n_boot block-wise bootstrapped samples.
#' @inheritParams bootstrap
#' @param balanced A boolean indicating if a balanced bootstrap procedure is
#' performed or not (default is TRUE).
#' @param keep_all_variables A boolean indicating if all variables have to be
#' kept even when some of them have null variance for at least one bootstrap
#' sample (default is FALSE).
#' @param pval For all the variables, a threshold for the proportion of the most
#' frequent value of this variable is computed. This threshold is evaluated so
#' that the probability to sample only this value is below `pval`. This
#' probability is corrected by the number of bootstrap samples and the number
#' of variables (default is `1e-15`).
#' @return \item{full_idx}{A list of size n_boot containing the observations
#' kept for each bootstrap sample.}
#' @return \item{summarize_column_sd_null}{A list of size the number of block
#' containing the variables that were removed from each block for all the
#' bootstrap samples. Variables are removed if they appear to be of null variance
#' in at least one bootstrap sample. If no variable is removed, return NULL.}
#' @title Generate sampled blocks for bootstrap.
#' @keywords internal

generate_resampling <- function(rgcca_res, n_boot, balanced = TRUE,
                                keep_all_variables = FALSE, pval = 1e-15){
  # Initialization
  NO_null_sd_var      = FALSE
  iter                = 0
  raw_blocks          = rgcca_res$call$raw
  N                   = NROW(raw_blocks[[1]])
  prob                = rep(1/N, N)
  # For any variable, threshold for the proportion of the most frequent value
  # of a variable. This threshold is computed so that the probability to sample
  # only this value is below `pval`. This probability is corrected by the number
  # of bootstrap samples and the number of variables.
  risky_threshold = max(1/N, (pval/(n_boot*sum(sapply(raw_blocks, NCOL))))^(1/N))
  # Identify variables with value having an observed proportion higher than
  # risky_threshold.
  risky_var = lapply(raw_blocks,
                     function(block)
                       which(apply(block, 2,
                                   function(x) max(table(x)/N) > risky_threshold)))
  # Keep only risky variables for each block.
  raw_blocks_filtered = mapply(function(x, y) x[, y, drop = FALSE],
                               raw_blocks, risky_var)
  # While there are variables with null variance among the risky variables.
  while (!NO_null_sd_var){
    if (balanced){# Balanced bootstrap sampling.
      full_idx = rep(x = seq(N), each = n_boot)
      full_idx = sample(x = full_idx, size = N*n_boot, replace = FALSE)
      full_idx = split(full_idx, ceiling(seq_along(full_idx)/N))
    }else{# Unbalanced bootstrap sampling.
      full_idx = lapply(seq(n_boot), function(x) sample(x = seq(N), replace = TRUE, prob = prob))
    }
    # Compute blocks for each bootstrap sample (only with risky variables).
    boot_blocks_filtered = lapply(full_idx, function(idx) lapply(raw_blocks_filtered,
                                                                 function(x){
                                                                   y           = x[idx, , drop = FALSE]
                                                                   rownames(y) = paste("S",1:length(idx))
                                                                   return(y)}))
    # For each sample, identify variables with a single unique value.
    boot_column_sd_null  = lapply(boot_blocks_filtered,
                                  function(boot)
                                    lapply(boot, function(boot_bl)
                                      which(apply(boot_bl, 2, function(x) length(unique(x)) == 1))))
    # Summarize through all the samples.
    eval_boot_sample    = sapply(boot_column_sd_null, function(x) sum(sapply(x, length)))
    NO_null_sd_var      = (sum(eval_boot_sample) == 0)
    if (NO_null_sd_var){ # If through all samples, all variables have non null variances.
      summarize_column_sd_null = NULL
    }else{# If at least one sample have been identified with a null variance variable.
      if (!keep_all_variables){# It IS allowed to remove variables.
        # Extract the troublesome variables.
        summarize_column_sd_null = Reduce("rbind", boot_column_sd_null)
        rownames(summarize_column_sd_null) = NULL
        summarize_column_sd_null = apply(summarize_column_sd_null, 2,
                                         function(x) unique(names(Reduce("c", x))))
        summarize_column_sd_null = mapply(function(x, y){
          z        = match(x, y)
          names(z) = x
          return(z)
        }, summarize_column_sd_null, lapply(raw_blocks, colnames))
        # Check if a whole block is troublesome
        is_full_block_removed = mapply(function(x, y) dim(x)[2] == length(y),
                                       raw_blocks,
                                       summarize_column_sd_null,
                                       SIMPLIFY = "array")
        if (sum(is_full_block_removed) == 0){# A whole block is NOT troublesome
          NO_null_sd_var = TRUE
          warning(paste("Variables: ",
                        paste(names(Reduce("c", summarize_column_sd_null)), collapse = " - "),
                        "appear to be of null variance in some bootstrap samples",
                        "and thus were removed from all samples. \n",
                        "==> RGCCA is run again without these variables."))
        }else{# A whole block IS troublesome
          # Generate at most five different re-sampling until is not the case anymore.
          if (iter > 5){# Otherwise STOP.
            stop_rgcca(paste("The variance of all the variables from blocks: ",
                             paste(names(raw_blocks)[is_full_block_removed], collapse = " - "),
                             "appear to be null in some bootstrap samples.",
                             "Please consider removing them."))
          }else{
            iter = iter + 1
          }
        }
      }else{# It is NOT allowed to remove variables.
        # Generate at most five different re-sampling until not a single variable
        # has a null variance.
        if (iter > 5){# Otherwise STOP.
          # Extract the troublesome variables.
          summarize_column_sd_null = Reduce("rbind", boot_column_sd_null)
          rownames(summarize_column_sd_null) = NULL
          summarize_column_sd_null = apply(summarize_column_sd_null, 2,
                                           function(x) unique(names(Reduce("c", x))))
          summarize_column_sd_null = mapply(function(x, y){
            z        = match(x, y)
            names(z) = x
            return(z)
          }, summarize_column_sd_null, lapply(raw_blocks, colnames))
          error_message = paste("Impossible to define all bootstrap samples",
                                "without variables with null variance. Please",
                                "consider removing these variables: ",
                                paste(names(Reduce("c", summarize_column_sd_null)), collapse = " - "))
          # In the balanced case, you CANNOT play with the sampling probability
          # of the different observations as it is unbalanced.
          if (balanced){
            error_message = paste0(error_message,
                                   ". Please, consider unbalanced bootstrap by",
                                   " setting 'balanced' to FALSE.")
          }
          stop_rgcca(error_message)
        }
        # In the unbalanced case, you CAN play with the sampling probability
        # of the different observations.
        if (!balanced){
          if (iter == 0){# The first time, you define your unbalancedness.
            # Each observed value of the risky variables is replaced by
            # `1 - the proportion of this observed value`, normalized so that
            # the sum through all the observations (for each variable) equals `1`.
            prob = sapply(raw_blocks_filtered,
                          function(block) apply(block, 2, function(var) {
                            occurences = table(var, useNA = "ifany")/length(var)
                            new_idx    = match(as.character(var), names(occurences))
                            new_var    = as.matrix(occurences[new_idx])
                            new_var    = (1-new_var)/sum(1-new_var)
                            return(new_var)
                          }))
            # The sampling probability for each observation is associated with
            # the maximum value of the previous matrix through all risky
            # variables (again normalize so that `sum(prob) = 1`). Thus
            # the sampling probability of an observation is more of less associated
            # with `1 - its proportion in the variable where this observation is in
            # the lowest frequent group (through all risky variables)`.
            prob = apply(Reduce('cbind', prob), 1, max)/sum(apply(Reduce('cbind', prob), 1, max))
          }
        }
        iter = iter + 1
      }
    }
  }
  return(list(full_idx = full_idx, summarize_column_sd_null = summarize_column_sd_null))
}
