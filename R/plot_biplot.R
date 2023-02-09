#' Plot the samples in the component plane along with the variables
#' composing the components.
#'
#' @inheritParams plot.rgcca
#' @param df Data frame with the data to plot.
#' @param theme_RGCCA Theme of the plot.
#' @noRd
plot_biplot <- function(df, title, x, block, comp, theme_RGCCA,
                        cex_point, sample_colors, sample_shapes,
                        show_labels, repel, var_colors, var_shapes,
                        show_arrows, ...) {
  sample_colors <- sample_colors[seq_along(levels(df$Y$response))]
  var_colors <- var_colors[seq_along(levels(df$a$response))]
  sample_shapes <- sample_shapes[seq_along(levels(df$Y$response))]
  var_shapes <- var_shapes[seq_along(levels(df$a$response))]

  names(sample_colors) <- levels(df$Y$response)
  names(var_colors) <- levels(df$a$response)
  names(sample_shapes) <- levels(df$Y$response)
  names(var_shapes) <- levels(df$a$response)

  # Prepare guide to show legend only if there is more than one block
  if (length(var_colors) > 1) {
    guide <- ggplot2::guide_legend(override.aes = list(
      colour = var_colors, shape = var_shapes
    ), order = 2)
  } else {
    guide <- "none"
    var_colors[1] <- "black"
  }

  colors <- c(sample_colors, var_colors)
  shapes <- c(sample_shapes, var_shapes)
  # Construct sample plot
  p <- plot_sample(df$Y, title, x, block, comp, theme_RGCCA,
                   cex_point, colors, shapes,
                   show_labels, repel)

  if (show_labels) {
    if (repel) {
      p <- p + geom_text_repel(
        data = df$a, aes(
          label = rownames(df$a), x = df$a[, 1],
          y = df$a[, 2], alpha = .data$response,
          color = .data$response
        ),
        size = cex_point, show.legend = FALSE, hjust = 0.5, vjust = -1
      )
    } else {
      p <- p + ggplot2::geom_text(
        data = df$a, aes(
          label = rownames(df$a), x = df$a[, 1],
          y = df$a[, 2], alpha = .data$response,
          color = .data$response
        ),
        size = cex_point, show.legend = FALSE, hjust = 0.5, vjust = -1
      )
    }
  }

  p <- p +
    ggplot2::geom_point(data = df$a, aes(
      x = df$a[, 1], y = df$a[, 2], alpha = .data$response,
      shape = .data$response, color = .data$response
    ), size = .5 * cex_point)

  if (show_arrows) {
    p <- p +
      ggplot2::geom_segment(
        data = df$a, aes(
          x = 0, y = 0, xend = df$a[, 1], alpha = .data$response,
          yend = df$a[, 2], color = .data$response
        ), show.legend = FALSE,
        arrow = ggplot2::arrow(length = ggplot2::unit(.7 * cex_point, "mm"))
      )
  }

  # Create second legend
  p <- p +
    ggplot2::scale_alpha_manual(
      "Block", values = rep(1, length(var_colors)),
      guide = guide
    )

  return(p)
}
