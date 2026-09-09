#' @export
#' @rdname plot
#' @order 4
plot.rgcca_bootstrap <- function(x, block = seq_along(x$rgcca$call$blocks),
                                 comp = 1, type = c("weights", "loadings","factors"),
                                 empirical = TRUE, n_mark = 30,
                                 display_order = TRUE,
                                 show_stars = TRUE, title = NULL,
                                 cex = 1, cex_sub = 12 * cex,
                                 cex_main = 14 * cex, cex_lab = 12 * cex,
                                 cex_point = 3 * cex, colors = NULL,
                                 adj.method = "fdr", ...) {
#  ### Perform checks and parse arguments
#  stopifnot(is(x, "rgcca_bootstrap"))
#  type <- type[1]
#  type <- match.arg(type, c("weights", "loadings"))
#  lapply(block, function(i) check_blockx("block", i, x$rgcca$call$blocks))
#  Map(function(y, z) check_compx(y, y, x$rgcca$call$ncomp, z), comp, block)
#  check_integer("n_mark", n_mark)
#  colors <- check_colors(colors, type = "variables")
#
#  ### Build data frame
#  column_names <- columns <- c(
#    "estimate", "mean", "sd", "lower_bound", "upper_bound", "pval"
#  )
#  if (!empirical) {
#    columns <- c(
#      "estimate", "mean", "sd", "th_lower_bound", "th_upper_bound", "th_pval"
#    )
#  }
#  df <- x$stats[x$stats$type == type, ]
#  col_pval <- ifelse(empirical, "pval", "th_pval")
#  df[, col_pval] <- p.adjust(df[, col_pval], method = adj.method)
#  df <- df[df$block %in% names(x$rgcca$blocks)[block], ]
#  df <- df[df$comp == comp, ]
#  rownames(df) <- df$var
#
#  df <- df[, columns]
#
#  colnames(df) <- column_names
#
#  df <- df[unlist(lapply(x$rgcca$blocks[block], colnames)), ]
#
#  df$response <- as.factor(unlist(lapply(block, function(j) {
#    rep(names(x$rgcca$blocks)[j], NCOL(x$rgcca$blocks[[j]]))
#  })))
#
#  if (display_order) {
#    df <- df[order(abs(df$estimate), decreasing = TRUE), ]
#  }
# # df <- df[df[, "sd"] != 0, ]
#  
#  n_mark <- min(n_mark, NROW(df))
#  df <- data.frame(df, order = seq(NROW(df), 1))[seq(n_mark), ]
#
#
# 
#
#  significance <- rep("", n_mark)
#  significance[df$pval < 1e-3] <- "***"
#  significance[df$pval >= 1e-3 & df$pval < 1e-2] <- "**"
#  significance[df$pval >= 1e-2 & df$pval < 5e-2] <- "*"
#
#  df$sign <- factor(significance,
#    levels = c(labels = c("***", "**", "*", ""))
#  )
#
#  ### Prepare plot
#  block_name <- ifelse(
#    length(block) > 1,
#    "",
#    paste0("(", names(x$rgcca$blocks)[block], ")")
#  )
#
#  title <- ifelse(is.null(title),
#    paste0(
#      "Bootstrap confidence interval ",
#      block_name, "\n (", type, " - ", x$n_boot,
#      " bootstrap samples - comp ", comp, ")"
#    ),
#    title
#  )
#
#  # Duplicate colors to avoid insufficient values in manual scale
#  colors <- rep(colors, NROW(df) / length(colors) + 1)
#
#  ### Construct plot
#  if (length(block) > 1) {
#    p <- ggplot(
#      df,
#      aes(x = .data$order, y = .data$estimate, color = .data$response)
#    ) +
#      ggplot2::scale_color_manual(values = colors) +
#      ggplot2::labs(color = "Block")
#  } else {
#    p <- ggplot(df, aes(x = .data$order, y = .data$estimate))
#  }
#
#  p <- p +
#    ggplot2::geom_point(size = .6 * cex_point) +
#    ggplot2::geom_errorbar(aes(
#      ymin = .data$lower_bound,
#      ymax = .data$upper_bound,
#      width = .2 * cex_point
#    ), linewidth = .2 * cex_point) +
#    ggplot2::coord_flip() +
#    theme_perso(cex, cex_main, cex_sub, cex_lab) +
#    ggplot2::labs(title = title, x = "", y = "") +
#    ggplot2::theme(
#      axis.line = ggplot2::element_blank(),
#      axis.ticks = ggplot2::element_blank(),
#      plot.margin = ggplot2::margin(5, 0, 0, 0, "mm"),
#      panel.grid.major.y = ggplot2::element_blank(),
#      panel.grid.minor.y = ggplot2::element_blank()
#    ) +
#    ggplot2::scale_x_continuous(
#      breaks = df$order,
#      labels = rownames(df)
#    ) +
#    ggplot2::geom_hline(
#      yintercept = 0, lty = "longdash", linewidth = .12 * cex_point
#    )
#  if (show_stars) {
#    p <- p + ggplot2::geom_text(
#      aes(label = .data$sign),
#      nudge_x = 0.1, size = 2 * cex_point, show.legend = FALSE
#    )
#  }
#
#  plot(p, ...)
#  invisible(p)
  stopifnot(is(x, "rgcca_bootstrap"))
  type <- type[1]
  type <- match.arg(type, c("weights", "loadings", "factors"))
  
  lapply(block, function(i) check_blockx("block", i, x$rgcca$call$blocks))
  Map(function(y, z) check_compx(y, y, x$rgcca$call$ncomp, z), comp, block)
  check_integer("n_mark", n_mark)
  colors <- check_colors(colors, type = "variables")

  ### Build data frame
  column_names <- columns <- c(
    "estimate", "mean", "sd", "lower_bound", "upper_bound", "pval"
  )
  if (!empirical) {
    columns <- c(
      "estimate", "mean", "sd", "th_lower_bound", "th_upper_bound", "th_pval"
    )
  }
  
  multi_blocks <- sapply(x$rgcca$call$blocks, function(b) {
    is.array(b) && length(dim(b)) > 2
  })
  
  blocks_names <- names(x$rgcca$blocks)
  
  df_list <- lapply(seq_along(block), function(ii) {
    j <- block[ii]  
    bn <- blocks_names[j]
    
    if ((length(dim(x$rgcca$blocks[[j]]))>2)& (type=='weights')){
      use_type='factors'
    
    }else{
      use_type=type
    }
    comp_j <- if (length(comp) == 1) comp else comp[ii]
    subset(x$stats, type == use_type & block == bn & comp == comp_j)
  })
  
  df <- do.call(rbind, df_list)
  if (is.null(df) || nrow(df) == 0) {
    message("No data to plot: no rows in bootstrap stats matching chosen blocks/type/comp.")
    return(invisible(NULL))
  }

  blkvec <- df$block
  
  col_pval <- ifelse(empirical, "pval", "th_pval")
  if (col_pval %in% colnames(df)) {
    df[, col_pval] <- p.adjust(df[, col_pval], method = adj.method)
  } else {
    df[, col_pval] <- NA
  }
  
  if (anyDuplicated(df$var) > 0) {
  # Fallback: df$var has duplicates, combine var and block
  rownames(df) <- paste(df$var, df$block, sep = "-")
} else {
  # Primary: df$var is completely unique
  rownames(df) <- df$var
}
  df <- df[, columns, drop = FALSE]
  colnames(df) <- column_names
  
  df$response <- factor(blkvec)
  if (sum(multi_blocks)>0){
  
    
  first=TRUE
    for (j in 1:length(block)){
      for (i in 1:length(dimnames(x$rgcca$blocks[[j]])[2:length( dimnames(x$rgcca$blocks[[j]]))]) ){
        if (first){
          y = rep(as.character(i),length(dimnames(x$rgcca$blocks[[j]])[2:length( dimnames(x$rgcca$blocks[[j]]))][i][[1]]) )
        first=FALSE
        }else{
          y = c(y,rep(as.character(i),length(dimnames(x$rgcca$blocks[[j]])[2:length( dimnames(x$rgcca$blocks[[j]]))][i][[1]]) ))

        }
      }
    }


    df$response=as.factor(y)


  }
  
  ### st

  if (display_order) {
    df <- df[order(abs(df$estimate), decreasing = TRUE), ]
  }
  df <- df[df[, "sd"] != 0, ]

  n_mark <- min(n_mark, NROW(df))
  df <- data.frame(df, order = seq(NROW(df), 1))[seq(n_mark), ]

  significance <- rep("", n_mark)
  significance[df$pval < 1e-3] <- "***"
  significance[df$pval >= 1e-3 & df$pval < 1e-2] <- "**"
  significance[df$pval >= 1e-2 & df$pval < 5e-2] <- "*"

  df$sign <- factor(significance,
    levels = c(labels = c("***", "**", "*", ""))
  )

  ### Prepare plot
  block_name <- ifelse(
    length(block) > 1,
    "",
    paste0("(", names(x$rgcca$blocks)[block], ")")
  )

  title <- ifelse(is.null(title),
    paste0(
      "Bootstrap confidence interval ",
      block_name, "\n (", type, " - ", x$n_boot,
      " bootstrap samples - comp ", comp, ")"
    ),
    title
  )

  # Duplicate colors to avoid insufficient values in manual scale
  colors <- rep(colors, NROW(df) / length(colors) + 1)
  
  ### Construct plot
  if (length(block) > 1) {
    p <- ggplot(
      df,
      aes(x = .data$order, y = .data$estimate, color = .data$response)
    ) +
      ggplot2::scale_color_manual(values = colors) +
      ggplot2::labs(color = "Block")
  } 
  else if (length(dim(x$rgcca$blocks[[block]]))>2){
        p <- ggplot(
      df,
      aes(x = .data$order, y = .data$estimate, color = .data$response)
    ) +
      ggplot2::scale_color_manual(values = colors) +
      ggplot2::labs(color = "Modes")

    }else {
    p <- ggplot(df, aes(x = .data$order, y = .data$estimate))
  }

  p <- p +
    ggplot2::geom_point(size = .6 * cex_point) +
    ggplot2::geom_errorbar(aes(
      ymin = .data$lower_bound,
      ymax = .data$upper_bound,
      width = .2 * cex_point
    ), linewidth = .2 * cex_point) +
    ggplot2::coord_flip() +
    theme_perso(cex, cex_main, cex_sub, cex_lab) +
    ggplot2::labs(title = title, x = "", y = "") +
    ggplot2::theme(
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(5, 0, 0, 0, "mm"),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor.y = ggplot2::element_blank()
    ) +
    ggplot2::scale_x_continuous(
      breaks = df$order,
      labels = rownames(df)
    ) +
    ggplot2::geom_hline(
      yintercept = 0, lty = "longdash", linewidth = .12 * cex_point
    )
  if (show_stars) {
    p <- p + ggplot2::geom_text(
      aes(label = .data$sign),
      nudge_x = 0.1, size = 2 * cex_point, show.legend = FALSE
    )
  }

  plot(p, ...)
  invisible(p)
}
