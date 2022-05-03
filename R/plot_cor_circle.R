#' Plot the block variables in the correlation circle
#'
#' @inheritParams plot.rgcca
#' @param df Data frame with the data to plot.
#' @param theme_RGCCA Theme of the plot.
#' @noRd
plot_cor_circle <- function(df, title, x, block, comp, theme_RGCCA,
                            cex_sub, cex_point, colors, shapes, ...) {
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
  title <- paste(title, names(x$call$blocks)[block[1]], sep = ": ")

  # Construct plot
  p <- ggplot(df, aes_(df[, 1], df[, 2], color = quote(response))) +
    geom_text(
      aes(label = rownames(df)),
      size = cex_point,
      show.legend = FALSE, hjust = 0.5, vjust = -1
    ) +
    theme_RGCCA +
    geom_vline(xintercept = 0, col = "grey", linetype = "dashed", size = 0.5) +
    geom_hline(yintercept = 0, col = "grey", linetype = "dashed", size = 0.5) +
    labs(title = title, x = xlab, y = ylab, color = "Block") +
    scale_y_continuous(breaks = NULL) +
    scale_x_continuous(breaks = NULL) +
    theme(
      legend.key.width = unit(nchar("Block"), "mm"),
      axis.text = element_blank(),
      axis.line = element_blank()
    )

  # Change colors and shapes based on discrete or continuous response
  discrete <- is.character(df$response) || is.factor(df$response)
  if (discrete) {
    p <- p + geom_point(aes_(shape = quote(response))) +
      scale_color_manual(values = colors) +
      scale_shape_manual(name = "Block", values = shapes)
  } else {
    p <- p + geom_point() +
      scale_color_gradient(low = colors[2], high = colors[3])
  }

  # Remove legend if response takes a single value
  if (length(unique(df$response)) == 1) {
    p <- p + theme(legend.position = "none")
  }

  p <- p +
    geom_path(
      aes_(quote(x), quote(y)),
      data = get_circle(),
      col = "grey",
      size =  0.5
    ) +
    geom_path(
      aes_(quote(x), quote(y)),
      data = get_circle() / 2,
      col = "grey",
      size = 0.5,
      lty = 2
    )

  return(p)
}
