#' Plot side by side the samples in the component plane and the
#' block variables in the correlation circle
#'
#' @inheritParams plot.rgcca
#' @param df Data frame with the data to plot.
#' @param theme_RGCCA Theme of the plot.
#' @noRd
plot_both <- function(df, title, x, block, comp, theme_RGCCA,
                      cex_point, sample_colors, sample_shapes,
                      show_sample_names, show_var_names,
                      repel, var_colors, var_shapes,
                      ...) {
  p <- grid.arrange(
    plot_sample(
      df[[1]], "Sample space", x, block, comp, theme_RGCCA,
      cex_point, sample_colors, sample_shapes, show_sample_names,
      repel, var_colors, var_shapes
    ),
    plot_cor_circle(
      df[[2]], "Correlation circle", x, block, comp, theme_RGCCA,
      cex_point, var_colors, var_shapes, show_var_names, repel
    ),
    nrow = 1, ncol = 2, top = title
  )

  invisible(p)
}
