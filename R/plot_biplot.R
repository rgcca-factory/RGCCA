#' Plot the samples in the component plane along with the variables
#' composing the components.
#'
#' @inheritParams plot.rgcca
#' @param df Data frame with the data to plot.
#' @param theme_RGCCA Theme of the plot.
#' @noRd
plot_biplot <- function(df, title, x, block, comp, theme_RGCCA,
                        cex_point, sample_colors, sample_shapes,
                        show_labels, repel, var_shapes, show_arrows,
                        ...) {
  # Construct sample plot
  p <- plot_sample(df$Y, title, x, block, comp, theme_RGCCA,
                   cex_point, sample_colors, sample_shapes,
                   show_labels, repel)

  df_a <- df$a

  if (show_labels) {
    if (repel) {
      p <- p + geom_text_repel(
        data = df$a,
        aes(label = rownames(df$a), x = df$a[, 1], y = df$a[, 2]),
        size = cex_point, color = "black",
        show.legend = FALSE, hjust = 0.5, vjust = -1
      )
    } else {
      p <- p + ggplot2::geom_text(
        data = df$a,
        aes(label = rownames(df$a), x = df$a[, 1], y = df$a[, 2]),
        size = cex_point, color = "black",
        show.legend = FALSE, hjust = 0.5, vjust = -1
      )
    }
  }

  p <- p +
    ggplot2::geom_point(data = df$a, aes(
      x = df$a[, 1], y = df$a[, 2]
    ), size = .5 * cex_point, shape = var_shapes[1], color = "black")

  if (show_arrows) {
    p <- p +
      ggplot2::geom_segment(
        data = df$a,
        aes(x = 0, y = 0, xend = df$a[, 1], yend = df$a[, 2]),
        arrow = ggplot2::arrow(length = ggplot2::unit(cex_point, "mm")),
        color = "black"
      )
  }

  return(p)
}
