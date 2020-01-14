#' Plot a bootstrap
#'
#' Plot the top variables from a bootstrap
#'
#' @inheritParams plot_histogram
#' @inheritParams plot_var_2D
#' @param show_boot A boolean to show the bootstrap mean and sd on the graphic
#' @examples
#' library(RGCCA)
#' data("Russett")
#' blocks = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
#'     politic = Russett[, 6:11] )
#' rgcca_out = rgcca.analyze(blocks)
#' boot = bootstrap(rgcca_out, 2, FALSE, n_cores = 1)
#' selected.var = get_bootstrap(rgcca_out, boot, n_cores = 1)
#' plot_bootstrap(selected.var, rgcca_out)
#' @export
plot_bootstrap <- function(
    df,
    rgcca,
    show_boot = TRUE,
    n_mark = 30,
    cex = 1,
    cex_sub = 16 * cex,
    cex_axis = 10 * cex) {

    color <- intneg <- intpos <- NULL
    J <- names(rgcca$a)

    if (NROW(df) > n_mark)
        df <- df[seq(n_mark), ]

    if (rgcca$superblock) {
        color2 <- factor(df$color)
        levels(color2) <- color_group(color2)
        p <- ggplot(df, aes(order, mean, fill = color))
    } else{
        color2 <- "black"
        p <- ggplot(df, aes(order, mean, fill = abs(mean)))
    }

    p <- plot_histogram(
        p,
        df,
        "Variable mean",
        as.character(color2),
        cex = cex,
        cex_sub = cex_sub,
        cex_axis = cex_axis)

    if (show_boot) {
        p <- p +
            geom_line(aes(x = order, y = mean), inherit.aes = FALSE, lwd = 0.7) +
            geom_point(aes(x = order, y = mean), inherit.aes = FALSE, size = 1.5)

        if (is(rgcca, "rgcca" ))
            p <- p +
                geom_errorbar(aes(ymin = intneg, ymax = intpos))
    }

    if (rgcca$superblock)
        col <- J
    else
        col <- J[-length(J)]

    if (rgcca$superblock) {
        matched <- match(rev(unique(df$color)), col)
        p <- order_color(rgcca$a, p, matched, rgcca$superblock)
    }

    return(p)
}
