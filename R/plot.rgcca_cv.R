#' @export
#' @rdname plot
#' @order 2
plot.rgcca_cv <- function(x, type = c("sd", "quantile"),
                          cex = 1,
                          cex_main = 14 * cex,
                          cex_sub = 12 * cex,
                          cex_point = 3 * cex,
                          cex_lab = 12 * cex,
                          display_order = TRUE, ...) {
  ### Perform checks and parse params
  stopifnot(is(x, "rgcca_cv"))
  type <- type[1]
  type <- match.arg(type, c("quantile", "sd"))

  ### Build data frame
  df <- data.frame(
    combinations = x$stats$combinations,
    ymin = apply(x$cv, 1, min),
    ymax = apply(x$cv, 1, max)
  )
  if (type == "quantile") {
    df$middle <- x$stats$median
    df$lower <- x$stats$Q1
    df$upper <- x$stats$Q3
  } else {
    df$middle <- x$stats$mean
    df$lower <- x$stats$mean - x$stats$sd
    df$upper <- x$stats$mean + x$stats$sd
  }

  best <- which(apply(
    x$params, 1, function(z) identical(z, drop(x$best_params))
  ))

  idx_order <- seq_len(nrow(df))
  if (display_order) {
    idx_order <- sort(df$middle, decreasing = FALSE, index.return = TRUE)$ix
    df <- df[idx_order, ]
    best <- which(idx_order == best)
  }

  df$category <- rep("param", NROW(x$cv))
  df$category[best] <- "best_param"

  labels <- as.expression(df$combinations)
  labels[[best]] <- bquote(underline(bold(.(labels[[best]]))))

  df$combinations <- factor(
    df$combinations, levels = df$combinations, ordered = TRUE
  )

  df_points <- data.frame(
    combinations = rep(df$combinations, NCOL(x$cv)),
    y = c(x$cv[idx_order, ]),
    category = rep(df$category, NCOL(x$cv))
  )

  ### Prepare plot
  n_run_str <- ifelse(
    x$n_run > 1, paste0(" and ", x$n_run, " run(s)"), ""
  )
  validation_str <- ifelse(
    x$validation == "kfold",
    paste0("kfold: with ", x$k, " folds", n_run_str, ""),
    "leave-one-out"
  )

  title <- paste0(
    "Cross-validated ", x$metric, " (", validation_str, ")\nBest parameters: ",
    df$combinations[best]
  )
  xlab <- paste0("Tuning parameter sets (", x$par_type, ")")
  ylab <- paste("Mean", x$metric)

  ### Construct plot
  p <- ggplot(data = df, aes(
    x = .data$combinations,
    group = .data$combinations,
    color = .data$category
  )) +
    ggplot2::coord_flip() +
    ggplot2::scale_color_manual(
      values = c("param" = "grey", "best_param" = "black")
    ) +
    ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab) +
    ggplot2::scale_x_discrete(
      labels = labels, breaks = df$combinations,
      guide = ggplot2::guide_axis(check.overlap = TRUE)
    ) +
    ggplot2::geom_boxplot(
      aes(
        ymin = .data$ymin, lower = .data$lower, middle = .data$middle,
        upper = .data$upper, ymax = .data$ymax
      ),
      stat = "identity"
    ) +
    ggplot2::geom_jitter(
      data = df_points, aes(
        x = .data$combinations, y = .data$y, color = .data$category
      ), size = .5 * cex_point,
      position = ggplot2::position_jitter(height = 0.01, width = 0.2)
    ) +
    ggplot2::ggtitle(title) +
    theme_perso(cex, cex_main, cex_sub, cex_lab) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = .8 * cex_sub),
      legend.position = "none"
    )
  plot(p, ...)
  invisible(p)
}
