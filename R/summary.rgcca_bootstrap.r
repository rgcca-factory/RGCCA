#' @export
#' @rdname summary
#' @order 4
summary.rgcca_bootstrap <- function(object,
                                    block = seq_along(object$rgcca$call$blocks),
                                    comp = 1, type = c("weights", "loadings","factors"),
                                    empirical = TRUE, display_order = FALSE,
                                    adj.method = "fdr", ...) {
  ### Perform checks and parse arguments
  stopifnot(is(object, "rgcca_bootstrap"))
  type <- type[1]
  type <- match.arg(type, c("weights", "loadings","factors"))
  lapply(block, function(i) check_blockx("block", i, object$rgcca$call$blocks))
  Map(function(y, z) check_compx(y, y, object$rgcca$call$ncomp, z), comp, block)

  ### Construct data.frame
  column_names <- columns <- c(
    "estimate", "mean", "sd", "lower_bound",
    "upper_bound", "bootstrap_ratio", "pval", "adjust.pval"
  )
  if (!empirical) {
    columns <- c(
      "estimate", "mean", "sd", "th_lower_bound",
      "th_upper_bound", "bootstrap_ratio", "th_pval", "adjust.pval"
    )
  }
 # df <- object$stats[object$stats$type == type, ]
  multi_blocks <- sapply(object$rgcca$call$blocks, function(b) {
    is.array(b) && length(dim(b)) > 2
  })
  
  blocks_names <- names(object$rgcca$blocks)
  
  df_list <- lapply(seq_along(block), function(ii) {
    j <- block[ii]  
    bn <- blocks_names[j]
    if ((length(dim(object$rgcca$call$blocks[[j]]))>2)& (type=='weights')){
      use_type='factors'
    
    }else{
      use_type=type
    }
    comp_j <- if (length(comp) == 1) comp else comp[ii]
    subset(object$stats, type == use_type & block == bn & comp == comp_j)
  })
  
  df <- do.call(rbind, df_list)
  col_pval <- ifelse(empirical, "pval", "th_pval")
  if (col_pval %in% colnames(df)) {
    df[, col_pval] <- p.adjust(df[, col_pval], method = adj.method)
  } else {
    df[, col_pval] <- NA
  }
  df["adjust.pval"] <- p.adjust(df[, col_pval], method = adj.method)
 # df <- df[df$block %in% names(object$rgcca$blocks)[block], ]
 # df <- df[df$comp == comp, ]
  rownames(df) <- df$var
  df <- df[, columns]
  colnames(df) <- column_names

  if (display_order) {
    df <- df[order(abs(df$estimate), decreasing = TRUE), ]

  }# else {
    #df <- df[unlist(lapply(object$rgcca$blocks[block], colnames)), ]
  #}

  df <- format(df, digits = 3)

  ### Print
  print_call(object$rgcca$call)
  cat("\n")
  type_str <- ifelse(type == "weights", "weight", "loading")
  cat(paste0(
    "Extracted statistics from ", object$n_boot, " bootstrap samples.\n",
    "Block-", type_str, " vectors for component ", comp, ":"
  ), "\n")
  
  print(df, quote = FALSE, ...)
  cat("\n")
}
