# Default theme for ggplot
theme_perso <- function(cex = 1, cex_main = 25 * cex, cex_sub = 16 * cex) {
    theme(
        legend.text = element_text(size = 13 * cex),
        legend.title = element_text(face = "bold.italic", size = cex_sub ),
        plot.title = element_text(
            size = cex_main,
            face = "bold",
            hjust = 0.5,
            margin = margin(0, 0, 20, 0)
        )
    )
}
