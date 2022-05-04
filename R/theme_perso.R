# Default theme for ggplot
theme_perso <- function(cex = 1, cex_main = 12 * cex, cex_sub = 10 * cex,
                        cex_lab = 10 * cex) {
  ggplot2::theme_classic() +
    ggplot2::theme(
      legend.text = ggplot2::element_text(size = cex_lab),
      legend.title = ggplot2::element_text(
        face = "bold.italic", size = cex_sub
      ),
      plot.title = ggplot2::element_text(
        size = cex_main,
        face = "bold",
        hjust = 0.5,
        margin = ggplot2::margin(0, 0, 0, 0)
      )
    )
}
