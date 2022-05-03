#' Histogram of the weights or the loadings
#'
#' Histogram of the weights or the loadings sorted in decreasing order.
#' @inheritParams plot.rgcca
#' @param df Data frame with the data to plot.
#' @param theme_RGCCA Theme of the plot.
#' @noRd
plot_loadings <- function(df, title, x, block, comp, theme_RGCCA,
                          cex_sub, cex_point, colors, shapes, ...) {
  # Add colors depending on looking at superblock or regular blocks
  is_superblock <- block[1] > length(x$call$raw)
  if (is_superblock) {
    p <- ggplot(df, aes_(x = quote(x), y = quote(y), fill = quote(response))) +
      scale_fill_manual(values = colors) +
      labs(fill = "Block")
  } else {
    p <- ggplot(df, aes_(x = quote(x), y = quote(y), fill = quote(x))) +
      scale_fill_gradient(low = colors[2], high = colors[3])
  }
  # Construct plot
  p <- p + geom_bar(stat = "identity") +
    theme_RGCCA +
    labs(title = title, x = "", y = "") +
    scale_y_discrete(limits = rev) +
    theme(
      axis.text.y = element_text(
        size = cex_sub,
        face = "italic",
        color = "gray40"
      ),
      axis.text.x = element_text(
        size = cex_sub,
        face = "italic",
        color = "gray40"
      ),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      plot.margin = margin(5, 0, 0, 0, "mm")
    )

  # Hide legend if not superblock
  if (!is_superblock) {
    p <- p + theme(legend.position = "none")
  }
  return(p)
}
