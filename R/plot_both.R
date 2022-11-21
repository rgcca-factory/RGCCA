#' Plot side by side the samples in the component plane and the
#' block variables in the correlation circle
#'
#' @inheritParams plot.rgcca
#' @param df Data frame with the data to plot.
#' @param theme_RGCCA Theme of the plot.
#' @noRd
plot_both <- function(df, title, x, block, comp, theme_RGCCA,
                      cex_main, cex_sub, cex_point, colors,
                      shapes, show_labels, repel) {
  p <- grid.arrange(
    plot_sample(
      df[[1]], "Sample space", x, block, comp, theme_RGCCA,
      cex_main, cex_sub, cex_point, colors, shapes, show_labels,
      repel
    ),
    plot_cor_circle(
      df[[2]], "Correlation circle", x, block, comp, theme_RGCCA,
      cex_main, cex_sub, cex_point, colors, shapes, show_labels,
      repel
    ),
    nrow = 1, ncol = 2
  ) +
    theme_RGCCA

  invisible(p)
}
