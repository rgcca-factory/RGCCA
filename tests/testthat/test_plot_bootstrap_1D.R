#' Plot a bootstrap in 1D
#'
#' Histogram of the best variables of an RGCCA bootstrap with, on the x-axis,
#' the number of non-zero occurrences (SGCCA) or the bootstrap-ratio
#' (mean/sd; RCCA). The bars are colored according to the average weight of
#' the boostrap  (according to an ascending gradient from red to blue)
#' @inheritParams plot_histogram
#' @param b A matrix of boostrap
#' @param x A character for the column used in the plot
#' @param y A character for the column to color the bars
#' @param n An integer giving the number maximum of top variables
#' @examples
#' library(RGCCA)
#' data("Russett")
#' blocks = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
#'     politic = Russett[, 6:11] )
#' rgcca_out = rgcca.analyze(blocks, tau = 0.75, type = "sgcca")
#' boot = bootstrap(rgcca_out, 2, n_cores = 1)
#' selected.var = get_bootstrap(rgcca_out, boot, n_cores = 1)
#' plot_bootstrap_1D(selected.var)
#' rgcca_out = rgcca.analyze(blocks)
#' boot = bootstrap(rgcca_out, 2, n_cores = 1)
#' selected.var = get_bootstrap(rgcca_out, boot, n_cores = 1)
#' plot_bootstrap_1D(selected.var)
#' @export
plot_bootstrap_1D <- function(
    b,
    x = "br",
    y = "occ",
    n = 50,
    cex = 1,
    cex_sub = 16 * cex,
    cex_axis = 10 * cex) {

    set_occ <- function(x) {
        match.arg(x, names(attributes(b)$indexes))
        if (x == "occ" && !x %in% colnames(b))
            return("sign")
        else
            return(x)
    }

    x <- set_occ(x)
    y <- set_occ(y)

    if (y == "sign") 
        color = seq(2)
    else
        color = 1

    b <- head(b, n)
    p <- ggplot(
        b,
        aes(x = order,
            y = b[, x],
            fill = b[, y]))

    plot_histogram(
        p,
        b,
        attributes(b)$indexes[[x]],
        color,
        low_col = color_group(seq(3))[1],
        mid_col = "white",
        high_col = color_group(seq(3))[3],
        cex = cex,
        cex_sub = cex_sub,
        cex_axis = cex_axis) +
    labs(fill = attributes(b)$indexes[[y]])
}
