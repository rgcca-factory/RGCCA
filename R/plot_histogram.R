#' Histogram settings
#'
#' Default font for a vertical barplot.
#'
#' @inheritParams plot2D
#' @param group A vector of character giving the group for the rows
#' @param low_col A character giving the color used for the lowest part of
#' the gradient
#' @param high_col A character giving the color used for the highest part of
#' the gradient
#' @param mid_col A character giving the color used for the middle part of
#' the gradient
#' @param cex_axis An integer for the size of the axis text
#' @param colors reoresenting a vector of colors
#' @examples
#' df = data.frame(x = runif(30), order = 30:1)
#' library("ggplot2")
#' p = ggplot(df, aes(order, x))
#' plot_histogram(p, df, "This is my title")
#' # Add colors per levels of a variable
#' df$color = rep(c(1,2,3), each=10)
#' p = ggplot(df, aes(order, x, fill = color))
#' plot_histogram(p, df, "Histogram", as.character(df$color))
#' @export
plot_histogram <- function(
    p,
    df,
    title = "Histogram",
    group = NA,
    colors = NULL,
    cex = 1,
    cex_main = 25 * cex,
    cex_sub = 16 * cex,
    cex_axis = 10 * cex
) {

    colors <- check_colors(colors)

    if (NROW(df) <= 10 || is(df, "d_ave")) {
        width <- NULL
        if (is(df, "d_ave"))
            cex_axis <- 12
    } else
        width <- 1

    if (NROW(df) < 3)
        mar <- 60
    else if (NROW(df) < 5)
        mar <- 30
    else
        mar <- 0

    axis <- function(margin){
        element_text(
            size = cex_axis,
            face = "italic",
            color = "gray40"
        )
    }

    p <- p + geom_bar(stat = "identity", width = width) +
        coord_flip() + labs(title = title,  x = "", y = "") +
        theme_classic() +
        theme_perso(cex, cex_main, cex_sub) +
        theme(
            axis.text.y = axis(),
            axis.text.x = axis(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            plot.subtitle = element_text(
                hjust = 0.5,
                size = cex_sub,
                face = "italic"
            ),
            plot.margin = margin(0, 0, mar, 0, "mm")
    )

    if  (!is(df, "d_ave")) {
            p <- p +
                scale_x_continuous(breaks = df$order, labels = rownames(df)) +
                labs(fill = "Blocks")
            if (length(group) == 1) {
                if (is.na(colors[2])) {
                    p <- p +
                        scale_fill_gradient(low = colors[1], high = colors[3]) +
                        theme(legend.position = "none")
                }else
                    p <- p +
                        scale_fill_gradient2(low = colors[1], high = colors[3], mid = colors[2])
            }
    }

    return(p)
}
