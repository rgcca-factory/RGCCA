#' Bar plots of the weights or the loadings
#'
#' Bar plots of the weights or the loadings sorted in decreasing order.
#' @inheritParams plot.rgcca
#' @param df Data frame with the data to plot.
#' @param theme_RGCCA Theme of the plot.
#' @noRd
plot_loadings <- function(df, title, x, block, comp, theme_RGCCA,
                          cex_sub, cex_point, colors, shapes) {
  # Add colors depending on looking at superblock or regular blocks
  is_superblock <- block[1] > length(x$call$raw)
  if (is_superblock) {
    p <- ggplot(df, aes(x = .data$x, y = .data$y, fill = .data$response)) +
      ggplot2::scale_fill_manual(values = colors) +
      ggplot2::labs(fill = "Block")
  } else {
    p <- ggplot(df, aes(x = .data$x, y = .data$y, fill = .data$x)) +
      ggplot2::scale_fill_gradient(low = colors[2], high = colors[3])
  }
  # Construct plot
  p <- p + ggplot2::geom_bar(stat = "identity") +
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
  if (!is_superblock) {
    p <- p + ggplot2::theme(legend.position = "none")
  }
  return(p)
}
