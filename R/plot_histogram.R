#' Histogram settings
#'
#' Default font for a vertical barplot.
#'
#' @inheritParams plot2D
#' @param group A vector of character giving the group for the rows
#' @param cex_axis An integer for the size of the axis text
#' @importFrom ggplot2 ggplot
plot_histogram <- function(
    p,
    df,
    title = "",
    group = NA,
    colors = NULL,
    cex = 1,
    cex_main = 14 * cex,
    cex_sub = 12 * cex,
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
       mar <- 50
   else if (NROW(df) <= 15)
       mar <- round(3 / NROW(df) * 50)
   else
        mar <- 0
    print(mar)
    if (NROW(df) > 50)
        cex_axis <- 7
    if (NROW(df) > 75)
        cex_axis <- 5

    axis <- function(margin){
        element_text(
            size = cex_axis,
            face = "italic",
            color = "gray40"
        )
    }

    p <- p + geom_bar(stat = "identity", width = width) +
        coord_flip() + 
        labs(title = title,  x = "", y = "") +
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
                colors <- c(color_group(seq(3))[3],  "gray", color_group(seq(3))[1])
            p <- p +
                scale_fill_gradientn(colors = colors, na.value = "black")
        } else  if ((is.character2(group[!is.na(group)]) ||
                            length(unique(group)) <= 5 )) {

            cols=color_group(group, colors)
           p <- p + scale_fill_manual(values = cols,limit=names(cols),drop=FALSE)
        }
    }

    return(p)
}
