#' Histogram of Average Variance Explained
#'
#' Histogram of the model quality (based on Average Variance Explained)
#' for each blocks and sorted in decreasing order.
#' @inheritParams plot.rgcca
#' @param df Data frame with the data to plot.
#' @param theme_RGCCA Theme of the plot.
#' @noRd
plot_ave <- function(df, title, x, block, comp, theme_RGCCA,
                     cex_sub, cex_point, colors, shapes, ...) {
  # Construct plot
  p <- ggplot(df, aes_(x = quote(AVE), y = quote(block), fill = quote(comp))) +
    geom_bar(position = position_stack(reverse = TRUE), stat = "identity") +
    stat_identity(
      geom = "text", colour = "black", size = cex_point,
      aes_(label = quote(AVE)),
      position = position_stack(reverse = TRUE, vjust = 0.5)
    ) +
    scale_fill_manual(values = colors) +
    theme_RGCCA +
    labs(subtitle = print_comp(x, outer = TRUE)) +
    labs(title = title, x = "", y = "", fill = "Component") +
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
      plot.subtitle = element_text(
        hjust = 0.5,
        size = cex_sub,
        face = "italic"
      ),
      plot.margin = margin(5, 0, 0, 0, "mm")
    )
  return(p)
}
