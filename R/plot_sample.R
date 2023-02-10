#' Plot the samples in the component plane
#'
#' @inheritParams plot.rgcca
#' @param df Data frame with the data to plot.
#' @param theme_RGCCA Theme of the plot.
#' @noRd
plot_sample <- function(df, title, x, block, comp, theme_RGCCA,
                        cex_point, sample_colors, sample_shapes,
                        show_sample_names, repel, var_colors,
                        var_shapes, ...) {
  xlab <- print_comp(x, comp[1], block[1])
  ylab <- print_comp(x, comp[2], block[2])

  discrete <- is.character(df$response) || is.factor(df$response)
  if (discrete) {
    df$response <- as.factor(df$response)
    sample_colors <- sample_colors[seq_along(levels(df$response))]
    sample_shapes <- sample_shapes[seq_along(levels(df$response))]
    names(sample_colors) <- levels(df$response)
    names(sample_shapes) <- levels(df$response)
  }

  # Change axis labels and title depending on using one are two blocks
  if (block[1] != block[2]) {
    xlab <- paste(xlab, names(x$blocks)[block[1]], sep = " - ")
    ylab <- paste(ylab, names(x$blocks)[block[2]], sep = " - ")
  } else {
    title <- paste(title, names(x$blocks)[block[1]], sep = ": ")
  }

  # Construct plot
  p <- ggplot(df, aes(df[, 1], df[, 2]))
  if (show_sample_names) {
    if (repel) {
      p <- p + geom_text_repel(
        aes(label = rownames(df), color = .data$response),
        size = cex_point,
        show.legend = FALSE, hjust = 0.5, vjust = -1
      )
    } else {
      p <- p + ggplot2::geom_text(
        aes(label = rownames(df), color = .data$response),
        size = cex_point,
        show.legend = FALSE, hjust = 0.5, vjust = -1
      )
    }
  }
  p <- p +
    theme_RGCCA +
    ggplot2::geom_vline(
      xintercept = 0, linetype = "dashed", linewidth = .1 * cex_point
    ) +
    ggplot2::geom_hline(
      yintercept = 0, linetype = "dashed", linewidth = .1 * cex_point
    ) +
    ggplot2::labs(title = title, x = xlab, y = ylab, color = "Response") +
    ggplot2::theme(
      legend.key.width = ggplot2::unit(nchar("Response"), "mm"),
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    ) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = .1)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = .1))

  # Change colors and shapes based on discrete or continuous response
  if (discrete) {
    # Remove legend if response takes a single value
    if (length(unique(df$response)) == 1) {
      guide <- "none"
    } else {
      guide <- ggplot2::guide_legend(order = 1)
    }
    p <- p +
      ggplot2::geom_point(
        aes(shape = .data$response, color = .data$response),
        size = .5 * cex_point
      ) +
      ggplot2::scale_color_manual(
        values = c(sample_colors, var_colors), guide = guide,
        breaks = df$response
      ) +
      ggplot2::scale_shape_manual(
        name = "Response", values = c(sample_shapes, var_shapes),
        guide = guide, breaks = df$response
      )
  } else {
    p <- p + ggplot2::geom_point(
      size = .5 * cex_point, aes(color = .data$response)
      ) +
      ggplot2::scale_color_gradient(
        low = sample_colors[1],
        high = sample_colors[2],
        guide = ggplot2::guide_colourbar(order = 1)
      )
  }

  return(p)
}
