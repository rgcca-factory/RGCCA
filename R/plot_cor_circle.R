#' Plot the block variables in the correlation circle
#'
#' @inheritParams plot.rgcca
#' @param df Data frame with the data to plot.
#' @param theme_RGCCA Theme of the plot.
#' @noRd
plot_cor_circle <- function(df, title, x, block, comp, theme_RGCCA,
                            cex_main, cex_sub, cex_point, colors,
                            shapes, show_labels, repel) {
  # Auxiliary function to construct circles
  get_circle <- function(center = c(0, 0), diameter = 2, npoints = 100) {
    r <- diameter / 2
    tt <- seq(0, 2 * pi, length.out = npoints)
    xx <- center[1] + r * cos(tt)
    yy <- center[2] + r * sin(tt)
    return(data.frame(x = xx, y = yy))
  }

  xlab <- print_comp(x, comp[1], block[1])
  ylab <- print_comp(x, comp[2], block[1])

  # Change axis labels and title depending on using one are two blocks
  title <- paste(title, names(x$blocks)[block[1]], sep = ": ")

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
    ggplot2::labs(title = title, x = xlab, y = ylab, color = "Block") +
    ggplot2::theme(
      legend.key.width = ggplot2::unit(nchar("Block"), "mm"),
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
      ggplot2::scale_color_manual(values = colors) +
      ggplot2::scale_shape_manual(name = "Block", values = shapes)
  } else {
    p <- p + ggplot2::geom_point(size = .5 * cex_point) +
      ggplot2::scale_color_gradient(low = colors[2], high = colors[3])
  }

  # Remove legend if response takes a single value
  if (length(unique(df$response)) == 1) {
    p <- p + ggplot2::theme(legend.position = "none")
  }

  p <- p +
    ggplot2::geom_path(
      aes(.data$x, .data$y),
      data = get_circle(),
      col = "black",
      linewidth = .1 * cex_point
    ) +
    ggplot2::geom_path(
      aes(.data$x, .data$y),
      data = get_circle() / 2,
      col = "black",
      linewidth = .1 * cex_point,
      linetype = "dashed"
    )

  return(p)
}
