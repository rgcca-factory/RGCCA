#' Plot a bootstrap in 2D
#'
#' Biplot of the top variables from a SGCCA bootstrap with the number of
#' non-zero occurences in x-axis and the boot-ratio (mean/sd) in y-axis.
#' Negative weights are colored in red and the positive ones are in green.
#'
#' @inheritParams plot2D
#' @inheritParams get_bootstrap
#' @param b A boostrap object \code{\link[RGCCA]{bootstrap}}
#' @param x A character for the index to plot in x-axis (see details).
#' @param y A character for the index to plot in y-axis (see details).
#' @param df_b A get_bootstrap object \code{\link[RGCCA]{get_bootstrap}}
#' @details 
#' \itemize{
#' \item 'estimate' for RGCCA weights
#' \item 'bootstrap_ratio' for the mean of the bootstrap weights / their standard error
#' \item 'sign' for significant 95% bootstrap interval
#' \item 'estimate' for RGCCA weights
#' \item 'occurrences' for non-zero occurences
#' \item 'mean' for the mean of the bootstrap weights
#' }
#' @examples
#' library(RGCCA)
#' data("Russett")
#' blocks = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
#'     politic = Russett[, 6:11] )
#' rgcca_out = rgcca(blocks, sparsity = 0.75, type = "sgcca")
#' boot = bootstrap(rgcca_out, 2, n_cores = 1)
#' plot_bootstrap_2D(boot, n_cores = 1)
#' rgcca_out = rgcca(blocks)
#' boot = bootstrap(rgcca_out, 2, n_cores = 1)
#' selected.var = get_bootstrap(boot, n_cores = 1,display_order=TRUE)
#' plot_bootstrap_2D(boot, n_cores = 1)
#' plot_bootstrap_2D(df_b = selected.var,n_cores=1)
#' @export
#' @seealso \code{\link[RGCCA]{bootstrap}}, \code{\link[RGCCA]{get_bootstrap}}
plot_bootstrap_2D <- function(
    b = NULL,
    df_b = NULL,
    x = "bootstrap_ratio",
    y = "occurrences",
    title = paste("Variable selection \nby",
           attributes(b)$n_boot,
           "bootstraps"),
    colors = NULL,
    cex = 1,
    cex_main = 14 * cex,
    cex_sub = 12 * cex,
    cex_point = 3 * cex,
    cex_lab = 10 * cex,
    comp = 1,
    i_block = NULL,
    collapse = FALSE,
    n_cores = parallel::detectCores() - 1) {

    if (missing(b) && missing(df_b))
        stop_rgcca("Please select a bootstrap object.")
    else if (!is.null(b)) {
        if (is.null(i_block))
            i_block <- length(b$bootstrap[[1]])
        df_b <- get_bootstrap(b, comp, i_block, collapse = collapse, n_cores = n_cores, display_order = TRUE)
    } else if (!is.null(df_b)) {
        if (is.null(i_block))
            i_block <- attributes(df_b)$n_blocks
        stopifnot(is(df_b, "df_bootstrap"))
    }

    title <- paste0(title, collapse = " ")
    check_ncol(list(df_b), 1)
    for (i in c("cex", "cex_main", "cex_sub", "cex_point", "cex_lab"))
        check_integer(i, get(i))
    check_colors(colors)

    set_occ <- function(x) {
        match.arg(x, names(attributes(df_b)$indexes))
        if (x == "occurrences" && !x %in% colnames(df_b))
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
            x[x == "NS"] <- 0
            x[x == "*"] <- 1
        }
        return(abs(as.double(x)))
    }

    p <- ggplot(
        df_b,
        aes(
            x = transform_x(df_b[, x]),
            y = transform_x(df_b[, y]),
            label = row.names(df_b),
            color = as.factor(!is.na(df_b[,"mean"]) & df_b[,"mean"] > 0  )
    )) +
    geom_text(size = cex_point * 0.75) +
    labs(
        y =  attributes(df_b)$indexes[[y]],
        x =  attributes(df_b)$indexes[[x]],
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
    scale_color_manual(values = as.vector(color_group(seq(2), colors = colors)))

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
