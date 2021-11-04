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
#' fit.stab = rgcca_stability(fit.sgcca, n_boot = 20, n_cores = 1)
#' boot = bootstrap(fit.stab,
#'                  n_boot = 30, n_cores = 1)
#' plot_bootstrap_1D(boot, type = "weight")
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
#' @importFrom grDevices grey.colors
plot_bootstrap_1D <- function(
    b = NULL,
    df_b = NULL,
    type = "weight",
    empirical = TRUE,
    x = "estimate",
    y = "sign",
    n_mark = 50,
    display_order = TRUE,
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
            empirical = empirical,
            display_order = display_order)
    }

    if (!is.null(df_b))
        stopifnot(is(df_b, "df_bootstrap"))

    check_integer("n_mark", n_mark)

    if (is.null(title)) {
        title <- paste0("Block-", type, " vector",
            "\n(", attributes(df_b)$n_boot,
            " bootstrap samples)")
    }

    lower_bound <- NULL -> upper_bound
    check_ncol(list(df_b), 1)

    if (n_mark > NROW(df_b)) n_mark <- NROW(df_b)

    if( display_order){
        df_b_head <- head(data.frame(
            order_df(df_b, "estimate", allCol = TRUE),
            order = NROW(df_b):1),
            n_mark)
    } else{
        df_b_head <- head(data.frame( df_b, order = NROW(df_b):1), n_mark)
    }

    df_b_head <- df_b_head[df_b_head[, "sd"] != 0, ]
    class(df_b_head) <- c(class(df_b), "d_boot1D")

    df_b_head$sign <- ""
    df_b_head$sign[df_b_head$pval<1e-3] <- "< 0.001"
    df_b_head$sign[df_b_head$pval>=1e-3 & df_b_head$pval<1e-2] <- "< 0.01"
    df_b_head$sign[df_b_head$pval>=1e-2 & df_b_head$pval<5e-2] <- "< 0.05"
    df_b_head$sign[df_b_head$pval>=5e-2 & df_b_head$pval<1e-1] <- "< 0.1"
    df_b_head$sign[df_b_head$pval>=1e-1] <- "> 0.1"

    df_b_head$sign = factor(df_b_head$sign,
                            levels=c(labels=c("< 0.001", "< 0.01",
                                              "< 0.05", "< 0.1", "> 0.1")))

    p <- ggplot(
        df_b_head,
        aes(x = order,
            y = df_b_head[, "estimate"],
            fill = df_b_head[, "sign"])
            )

    p <- plot_histogram(p, df_b_head, title) +
        scale_x_continuous(breaks = df_b_head$order,
                           labels = rownames(df_b_head)) +
        scale_fill_manual(values = grey.colors(6)[2:6],
                          labels = c("< 0.001", "< 0.01", "< 0.05",
                                     "< 0.1", "> 0.1"),
                          drop = FALSE,
                          name = "Signif.")


    if (nrow(df_b_head) <= 50){
        p <- p + geom_errorbar(aes(ymin = lower_bound,
                                   ymax = upper_bound,
                                   width = 0.5))
    }

    return(p)
}
