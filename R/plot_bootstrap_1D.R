#' Plot a fitted bootstrap object in 1D
#'
#' Display bootstrap confidence intervals for RGCCA or the number of non-zero
#' occurrences for SGCCA. The bars are colored according to the
#' significancy  block weight vectors ('*' or 'ns'; see 'pval' in details for
#' \code{\link[RGCCA]{get_bootstrap}}) for RGCCA and according to the
#' occurrences of non-zero weights for SGCCA. In SGCCA, the significant
#' variables are those above the three bars, respectively, with an alpha = 0.05
#' (dark red), 0.01 (red) and 0.001 (light red).
#' @inheritParams plot_histogram
#' @inheritParams get_bootstrap
#' @inheritParams plot_var_2D
#' @param df_b A get_bootstrap object \code{\link[RGCCA]{get_bootstrap}}
#' @param b A fitted bootstrap object \code{\link[RGCCA]{bootstrap}}
#' @param type Character string indicating the bootstrapped object to print:
#' block-weight vectors ("weight", default) or block-loading vectors
#' ("loadings").
#' @param x indicator used in the plot (see details).
#' @param y A character string indicating for the index to color the bars
#' (see details).
#' @param display_bar A logical value. If TRUE colorbar for significant
#' variables is displayed.
#' @param ... Other parameters (see plot_histogram)
#' @details
#' \itemize{
#' \item 'estimate' of the block weight vectors
#' \item 'bootstrap_ratio' of the block weight/loading vectors
#' \item 'sign' for significant 95% bootstrap interval
#' \item 'occurrences' number of for non-zero occurrences
#' \item 'mean'  mean of the bootstraped block weight vectors
#' }
#' @examples
#' data("Russett")
#' blocks = list(agriculture = Russett[, seq(3)],
#'               industry = Russett[, 4:5],
#'               politic = Russett[, 6:11])
#' fit.sgcca = rgcca(blocks, sparsity = 0.75, method = "sgcca")
#' fit.stab = rgcca_stability(fit.sgcca, n_cores = 1)
#' boot = bootstrap(fit.stab,
#'                  n_boot = 30, n_cores = 1)
#' plot_bootstrap_1D(boot, type = "loadings")
#' rgcca_out = rgcca(blocks)
#' boot = bootstrap(rgcca_out, 2, n_cores = 1)
#' selected.var = get_bootstrap(boot, display_order=TRUE)
#' plot_bootstrap_1D(boot)
#' plot_bootstrap_1D(df_b = selected.var)
#'
#' @export
#' @importFrom ggplot2 ggplot
#' @importFrom stats qbinom
#' @importFrom utils head
plot_bootstrap_1D <- function(
    b = NULL,
    df_b = NULL,
    type = "weight",
    x = "estimate",
    y = "sign",
    n_mark = 50,
    title = NULL,
    colors = NULL,
    comp = 1,
    display_bar = TRUE,
    i_block = length(b$bootstrap[[1]][[1]]),
    ...) {

    if (missing(b) && missing(df_b))
        stop_rgcca("Please select a fitted bootstrap object.")
    if (!is.null(b)) {
        df_b <- get_bootstrap(
            b,
            type = type,
            comp,
            block = i_block,
            display_order = TRUE
        )
    }

    if (!is.null(df_b))
        stopifnot(is(df_b, "df_bootstrap"))

    check_integer("n_mark", n_mark)

    if (is.null(title)) {
        title <- paste0(
            attributes(df_b)$indexes[[x]],
            "\n(",
            attributes(df_b)$n_boot,
            " bootstrap samples)")
    }

    if (is.null(colors)) {
        if (y != "sign")
            colors <- c(color_group(seq(3))[1], "gray", color_group(seq(3))[3])
        else
            colors <- c(color_group(seq(3))[1], color_group(seq(3))[3])
    }
    lower_bound <- NULL -> upper_bound
    check_ncol(list(df_b), 1)

    if (y == "sign")
        group <- c("NS","*")
    else
        group <- NA

    if (n_mark > NROW(df_b))
        n_mark <- NROW(df_b)

    df_b_head <- head(
        data.frame(
            order_df(df_b[, -NCOL(df_b)], x, allCol = TRUE),
            order = NROW(df_b):1),  n_mark)
    df_b_head <- df_b_head[df_b_head[, "sd"] != 0, ]
    class(df_b_head) <- c(class(df_b), "d_boot1D")

    if (!is.null(df_b_head$sign)) {
        df_b_head$sign[df_b_head$sign == 1] <- "*"
        df_b_head$sign[df_b_head$sign == 0] <- "ns"
    }

    p <- ggplot(
        df_b_head,
        aes(x = order,
            y = df_b_head[, x],
            fill = df_b_head[, y]))

    p <- plot_histogram(
        p,
        df_b_head,
        title,
        group,
        colors,
        ...) +
    labs(fill = attributes(df_b)$indexes[[y]])

    if (nrow(df_b_head) <= 50)
        p <- p +
            geom_errorbar(
                aes(
                    ymin = lower_bound,
                    ymax = upper_bound,
                    width = 0.5))

    return(p)
}
