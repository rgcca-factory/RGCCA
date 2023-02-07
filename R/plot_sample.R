#' Plot the samples in the component plane
#'
#' @inheritParams plot.rgcca
#' @param df Data frame with the data to plot.
#' @param theme_RGCCA Theme of the plot.
#' @noRd
plot_sample <- function(df, title, x, block, comp, theme_RGCCA,
                        cex_point, sample_colors, sample_shapes,
                        show_labels, repel, ...) {
  xlab <- print_comp(x, comp[1], block[1])
  ylab <- print_comp(x, comp[2], block[2])

  # Change axis labels and title depending on using one are two blocks
  if (block[1] != block[2]) {
    xlab <- paste(xlab, names(x$blocks)[block[1]], sep = " - ")
    ylab <- paste(ylab, names(x$blocks)[block[2]], sep = " - ")
  } else {
    title <- paste(title, names(x$blocks)[block[1]], sep = ": ")
  }

  # Construct plot
  p <- ggplot(df, aes(df[, 1], df[, 2], color = .data$response))
  if (show_labels) {
    if (repel) {
      p <- p + geom_text_repel(
        aes(label = rownames(df)),
        size = cex_point,
        show.legend = FALSE, hjust = 0.5, vjust = -1
      )
    } else {
      p <- p + ggplot2::geom_text(
        aes(label = rownames(df)),
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
  discrete <- is.character(df$response) || is.factor(df$response)
  if (discrete) {
    p <- p +
      ggplot2::geom_point(aes(shape = .data$response), size = .5 * cex_point) +
      ggplot2::scale_color_manual(values = sample_colors) +
      ggplot2::scale_shape_manual(name = "Response", values = sample_shapes)
  } else {
    p <- p + ggplot2::geom_point(size = .5 * cex_point) +
      ggplot2::scale_color_gradient(
        low = sample_colors[1], high = sample_colors[2]
      )
  }

  # Remove legend if response takes a single value
  if (length(unique(df$response)) == 1) {
    p <- p + ggplot2::theme(legend.position = "none")
  }

  return(p)
}
