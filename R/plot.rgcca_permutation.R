#' @export
#' @rdname plot
#' @order 3
plot.rgcca_permutation <- function(x,
                                   type = c("crit", "zstat"),
                                   cex = 1,
                                   title = NULL,
                                   cex_main = 14 * cex,
                                   cex_sub = 12 * cex,
                                   cex_point = 3 * cex,
                                   cex_lab = 12 * cex,
                                   display_order = TRUE,
                                   show_legend = FALSE, ...) {
  ### Perform checks and parse params
  stopifnot(is(x, "rgcca_permutation"))
  type <- type[1]
  match.arg(type, c("crit", "zstat"))

  ### Build data frame
  df <- data.frame(
    x = x$stats[, type],
    label = "Other parameter set",
    combinations = x$stats$combinations
  )

  # Reorder dataframe according to the quantity of interest
  idx_order <- seq_len(nrow(df))
  if (display_order) {
    idx_order <- sort(df$x, decreasing = FALSE, index.return = TRUE)$ix
    df <- df[idx_order, ]
  }

  # Mark the best parameter set
  best <- which(apply(
    x$params[idx_order, ], 1, function(z) identical(z, x$best_params)
  ))
  df$label[best] <- "Best parameter set"

  df$combinations <- factor(
    df$combinations,
    levels = rev(df$combinations), ordered = TRUE
  )

  ### Prepare plot
  crit_title <- ifelse(x$call$method %in% sparse_methods(),
    "SGCCA criterion",
    "RGCCA criterion"
  )
  xlab <- ifelse(type == "zstat", "Z-score", crit_title)
  ylab <- paste0("Tuning parameter sets (", x$par_type, ")")

  breaks <- rev(levels(df$combinations))
  labels <- as.expression(breaks)
  labels[[best]] <- bquote(underline(bold(.(labels[[best]]))))

  title <- ifelse(
    missing(title),
    paste0(
      "Permutation scores (", x$n_perms, " runs) \n Best parameters: ",
      df$combinations[best]
    ),
    title
  )

  ### Construct plot
  # Main plot (values of the quantity of interest per set of parameters)
  p <- ggplot(
    data = df, mapping = aes(x = .data$x, y = .data$combinations)
  ) +
    ggplot2::geom_point(
      aes(shape = .data$label, size = .data$label, color = .data$label)
    ) +
    ggplot2::labs(
      title = title,
      y     = ylab,
      x     = xlab
    ) +
    theme_perso(cex, cex_main, cex_sub, cex_lab) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = .8 * cex_sub)
    ) +
    ggplot2::scale_shape_manual("Non permuted",
      values = c(
        "Best parameter set" = 17,
        "Other parameter set" = 2
      ), guide = ggplot2::guide_legend(override.aes = list(
        color = c("black", "grey")
      ))
    ) +
    ggplot2::scale_size_manual("Non permuted",
      values = c(
        "Best parameter set" = 1.5 * cex_point,
        "Other parameter set" = .5 * cex_point
      )
    ) +
    ggplot2::scale_y_discrete(
      labels = labels, breaks = breaks,
      guide = ggplot2::guide_axis(check.overlap = TRUE)
    )

  # Second part of the plot (boxplots of the permuted criteria)
  if (type == "crit") {
    df2 <- data.frame(
      combinations = rep(df$combinations, NCOL(x$permcrit)),
      x = c(x$permcrit[idx_order, ]),
      label = rep(df$label, NCOL(x$permcrit))
    )
    p$layers <- c(
      ggplot2::geom_boxplot(
        data = df2,
        aes(x = .data$x, y = .data$combinations, color = .data$label),
        size = 0.8
      ),
      p$layers
    )
    p <- p + ggplot2::scale_color_manual("Permuted", values = c(
      "Best parameter set" = "black",
      "Other parameter set" = "grey"
    ))
  } else {
    p <- p + ggplot2::scale_color_manual("Non permuted", values = c(
      "Best parameter set" = "black",
      "Other parameter set" = "grey"
    ))
  }

  # Move legend
  if (!show_legend) {
    p <- p + ggplot2::theme(legend.position = "none")
  } else {
    p <- p +
      ggplot2::theme(
        legend.position = c(.9, 1.),
        legend.justification = c("right", "top")
      )
  }

  plot(p, ...)
  invisible(p)
}
