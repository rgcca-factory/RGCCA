#' Plot a bootstrap in 2D
#'
#' Biplot of the top variables from a SGCCA bootstrap with the number of
#' non-zero occurences in x-axis and the boot-ratio (mean/sd) in y-axis.
#' Negative weights are colored in red and the positive ones are in green.
#'
#' @inheritParams plot2D
#' @param colors reoresenting a vector of colors
#' @param b A matrix of boostrap
#' @param x A character for the column to plot in x-axis
#' @param y A character for the column to plot in y-axis
#' @examples
#' library(RGCCA)
#' data("Russett")
#' blocks = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
#'     politic = Russett[, 6:11] )
#' rgcca_out = rgcca(blocks, sparsity = 0.75, type = "sgcca")
#' boot = bootstrap(rgcca_out, 2, n_cores = 1)
#' selected.var = get_bootstrap(boot, n_cores = 1)
#' plot_bootstrap_2D(selected.var)
#' rgcca_out = rgcca(blocks)
#' boot = bootstrap(rgcca_out, 2, n_cores = 1)
#' selected.var = get_bootstrap(boot, n_cores = 1)
#' plot_bootstrap_2D(selected.var)
#' @export
#' @seealso \code{\link[RGCCA]{bootstrap}}, \code{\link[RGCCA]{get_bootstrap}}

plot_bootstrap_2D <- function(
    b,
    x = "bootstrap_ratio",
    y = "occurrences",
    title = paste("Variable selection \nby",
           attributes(b)$n_boot,
           "bootstraps"),
    colors = NULL,
    cex = 1,
    cex_main = 25 * cex,
    cex_sub = 16 * cex,
    cex_point = 3 * cex,
    cex_lab = 19 * cex){

    stopifnot(is(b, "df_bootstrap"))
    title <- paste0(title, collapse = " ")
    check_ncol(list(b), 1)
    for (i in c("cex", "cex_main", "cex_sub", "cex_point", "cex_lab"))
        check_integer(i, get(i))
    check_colors(colors)

    set_occ <- function(x) {
        match.arg(x, names(attributes(b)$indexes))
        if (x == "occurrences" && !x %in% colnames(b))
            return("sign")
        else
            return(x)
    }

    x <- set_occ(x)
    y <- set_occ(y)

    axis <- function(margin){
        element_text(
        face = "italic",
        size = cex_lab * 0.75,
        margin = margin
        )
    }

    transform_x <- function(x){
        if ("*" %in% x) {
            x[x == ""] <- 0
            x[x == "*"] <- 1
        }
        return(abs(as.double(x)))
    }

    p <- ggplot(
        b,
        aes(
            x = transform_x(b[, x]),
            y = transform_x(b[, y]),
            label = row.names(b),
            color = as.factor(mean > 0)
    )) +
    geom_text(size = cex_point * 0.75) +
    labs(
        y =  attributes(b)$indexes[[y]],
        x =  attributes(b)$indexes[[x]],
        title = title
    ) +
    theme_classic() +
    theme_perso(cex, cex_main, cex_sub) +
    theme(
        legend.position = "none",
        axis.title.y = axis(margin(0, 20, 0, 0)),
        axis.title.x = axis(margin(20, 0, 0, 0)),
        axis.text = element_text(size = 13 * cex)
    ) +
    scale_color_manual(values = color_group(seq(2), colors = colors))

    limites <- function(p, x){
        if (x %in% c("sign", "occurrences")) {
            axis <- deparse(substitute(x))
            func <- get(paste0(axis, "lim"))
            p <- p + func(0, 1)
            if (x == "sign") {
                p <- p + 
                    get(paste("scale", axis, "discrete", sep = "_"))(
                        labels = c("ns", "*"),
                        limits = c(0, 1)
                    )
            }
        }
        return(p)
    }

    p <- suppressMessages(limites(p, x))
    p <- suppressMessages(limites(p, y))

    return(p)
}
