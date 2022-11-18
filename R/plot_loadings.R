#' Bar plots of the weights or the loadings
#'
#' Bar plots of the weights or the loadings sorted in decreasing order.
#' @inheritParams plot.rgcca
#' @param df Data frame with the data to plot.
#' @param theme_RGCCA Theme of the plot.
#' @noRd
plot_loadings <- function(df, title, x, block, comp, theme_RGCCA,
                          cex_main, cex_sub, cex_point, colors,
                          shapes, show_labels) {
  # Add colors depending on looking at superblock or regular blocks
  is_multiblock <- (length(block) > 1) || (block == length(x$call$raw) + 1)
  if (is_multiblock) {
    p <- ggplot(df, aes(x = .data$x, y = .data$y, color = .data$response)) +
      ggplot2::scale_color_manual(values = colors) +
      ggplot2::labs(color = "Block")
  } else {
    p <- ggplot(df, aes(x = .data$x, y = .data$y))
  }
  # Construct plot
  p <- p +
    ggplot2::geom_point(size = .5 * cex_point) +
    ggplot2::geom_linerange(aes(
      xmin = 0,
      xmax = .data$x
    ), size = .2 * cex_point) +
    ggplot2::geom_vline(xintercept = 0, lty = "longdash") +
    theme_RGCCA +
    ggplot2::labs(title = title, x = "", y = "") +
    ggplot2::scale_y_discrete(limits = rev) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(
        size = cex_sub,
        face = "italic",
        color = "gray40"
      ),
      axis.text.x = ggplot2::element_text(
        size = cex_sub,
        face = "italic",
        color = "gray40"
      ),
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(5, 0, 0, 0, "mm")
    )

  # Hide legend if not superblock
  if (!is_multiblock) {
    p <- p + ggplot2::theme(legend.position = "none")
  }
  return(p)
}
