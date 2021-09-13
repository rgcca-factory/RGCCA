#' Plot a bootstrap in 2D
#'
#' Graph of the best variables from a bootstrap with, in x-axis, the number of
#' non-zero occurences (SGCCA) or the significant 95% bootstrap 
#' intervals (RGCCA; '*' or 'ns'; see 'p.vals' in details for 
#' \code{\link[RGCCA]{get_bootstrap}}). In in y-axis are the bootstrap-ratios (mean/sd) .
#' Negative weights are colored in red and the positive ones are in green.
#'
#' @inheritParams plot2D
#' @inheritParams plot_bootstrap_1D
#' @details 
#' \itemize{
#' \item 'estimate' for RGCCA weights
#' \item 'bootstrap_ratio' for the mean of the bootstrap weights / their standard error
#' \item 'sign' for significant 95% bootstrap interval
#' \item 'occurrences' for non-zero occurences
#' \item 'mean' for the mean of the bootstrap weights
#' }
#' @examples
#' data("Russett")
#' blocks = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
#'     politic = Russett[, 6:11] )
#' rgcca_out = rgcca(blocks, sparsity = 0.75, method = "sgcca")
#' boot = bootstrap(rgcca_stability(rgcca_out, n_cores = 1), 2, n_cores = 1)
#' plot_bootstrap_2D(boot)
#' rgcca_out = rgcca(blocks)
#' boot = bootstrap(rgcca_out, 2, n_cores = 1)
#' selected.var = get_bootstrap(boot,display_order=TRUE)
#' plot_bootstrap_2D(boot)
#' plot_bootstrap_2D(df_b = selected.var)
#' @export
#' @seealso \code{\link[RGCCA]{bootstrap}}, \code{\link[RGCCA]{get_bootstrap}}
plot_bootstrap_2D <- function(
    b = NULL,
    df_b = NULL,
    type = "weight",
    x = "bootstrap_ratio",
    y = "sign",
    title = NULL,
    colors = NULL,
    cex = 1,
    cex_main = 14 * cex,
    cex_sub = 12 * cex,
    cex_point = 3 * cex,
    cex_lab = 10 * cex,
    comp = 1,
    i_block = NULL) {

    if (missing(b) && missing(df_b))
        stop_rgcca("Please select a fitted bootstrap object.")
    else if (!is.null(b)) {
        if (is.null(i_block))
            i_block <- length(b$bootstrap$W[[1]])
        df_b <- get_bootstrap(
            b,
            type = type,
            comp,
            block = i_block,
            display_order = TRUE)
    } else if (!is.null(df_b)) {
        stopifnot(is(df_b, "df_bootstrap"))
    }

    if (is.null(title)) {
        title <- paste0("Block-", type, " vector: ", 
            attributes(df_b)$block,
            "\n(", attributes(df_b)$n_boot,
            " bootstrap samples)")
    }
    title <- paste0(title, collapse = " ")
    check_ncol(list(df_b), 1)
    for (i in c("cex_main", "cex_sub", "cex_point", "cex_lab"))
        check_integer(i, get(i))
    check_integer("cex", cex, float = TRUE)
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

    df_b <- df_b[df_b[, "sd"] != 0, ]
    df_b$sign <- df_b$pval
    df_b$sign[df_b$sign > 0.05] <- 1
    df_b$sign[df_b$sign <= 0.05] <- 0

    transform_x <- function(x){
        return(abs(as.double(x)))
    }

    p <- ggplot(
        df_b,
        aes(
            x = transform_x(df_b[, x]),
            y = transform_x(df_b[, y]),
            label = row.names(df_b),
            color = as.factor(df_b[, "sign"])
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
                        labels = c("*", "NS"),
                        limits = c(0, 1)
                    )
            }
        }
        return(p)
    }

    p <- suppressWarnings(suppressMessages(limites(p, x)))
    p <- suppressWarnings(suppressMessages(limites(p, y)))

    return(p)
}
