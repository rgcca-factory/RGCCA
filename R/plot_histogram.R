#' Histogram settings
#'
#' Default font for a vertical barplot.
#'
#' @inheritParams plot2D
#' @param group A vector of character giving the group for the rows
#' @param cex_axis An integer for the size of the axis text
#' @param colors reoresenting a vector of colors
#' @examples
#' df = data.frame(x = runif(30), order = 30:1)
#' library("ggplot2")
#' p = ggplot(df, aes(order, x))
#' plot_histogram(p, df, "This is my title")
#' # Add colors per levels of a variable
#' df$color = rep(c(1,2,3), each=10)
#' p = ggplot(df, aes(order, x, fill = as.character(color)))
#' plot_histogram(p, df, "Histogram", as.character(df$color))
#' @export
plot_histogram <- function(
    p,
    df,
    title = "",
    group = NA,
    colors = NULL,
    cex = 1,
    cex_main = 25 * cex,
    cex_sub = 16 * cex,
    cex_axis = 10 * cex
) {

    for (i in c("cex", "cex_main", "cex_sub", "cex_axis"))
        check_integer(i, get(i))

    stopifnot(is(p, "ggplot"))
    check_colors(colors)
    title <- paste0(title, collapse = " ")
    group <- as.vector(group)

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
        if (length(group) == 1){
            if (is.null(colors))
                colors <- c("blue", "gray", "#cd5b45")
            p <- p +
                scale_fill_gradientn(colors = colors, na.value = "black")
        } else  if ((is.character2(group[!is.na(group)]) ||
                            length(unique(group)) <= 5 )) {
            p <- p + scale_fill_manual(values = color_group(group, colors))
        }
    }

    return(p)
}
