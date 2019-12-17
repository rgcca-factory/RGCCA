#' Histogram of a fingerprint
#'
#' Histogram of the higher outer weight vectors for a component of a block 
#' (by default, the superblock or the last one) analysed by R/SGCCA
#'
#' @inheritParams plot_var_2D
#' @inheritParams plot_histogram
#' @param comp An integer giving the index of the analysis components
#' @param type A string giving the criterion to selects biomarkers : either 
#' "cor" for correlation between the component and the block
#' or "weight" for the weight of the RGCCA
#' @seealso \code{\link[RGCCA]{rgcca}}, \code{\link[RGCCA]{sgcca}}
#' @examples
#' weights = lapply(seq(3), function(x) matrix(runif(7*2), 7, 2))
#' for (i in seq(3))
#' row.names(weights[[i]]) <- paste0(letters[i],
#'      letters[seq(NROW(weights[[i]]))])
#' weights[[4]] = Reduce(rbind, weights)
#' rgcca_out = list(a = weights)
#' names(rgcca_out$a) = LETTERS[seq(4)]
#' rgcca_out$blocks = lapply(rgcca_out$a, t)
#' rgcca_out$superblock = TRUE
#' # With the 1rst component of the superblock
#' plot_var_1D(rgcca_out, 1, type = "weight")
#' # With the 2nd component of the 1rst block by selecting the ten higher weights
#' plot_var_1D(rgcca_out, 2, 10, 1, type = "weight")
#' library(RGCCA)
#' data("Russett")
#' blocks = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
#'     politic = Russett[, 6:11] )
#' rgcca_out = rgcca.analyze(blocks)
#' plot_var_1D(rgcca_out, collapse = TRUE)
#' @export
plot_var_1D <- function(
    rgcca,
    comp = 1,
    n_mark = 100,
    i_block = length(rgcca$a),
    type = "cor",
    collapse = FALSE,
    cex = 1,
    cex_sub = 16 * cex,
    cex_axis = 10 * cex) {

    df <- get_ctr2(
        rgcca = rgcca,
        compx = comp,
        compy = comp,
        i_block = i_block,
        type = type,
        n_mark = n_mark,
        collapse = collapse,
        remove_var = FALSE
    )

    if (i_block < length(rgcca$a) || is(rgcca, "pca"))
        rgcca$superblock <- FALSE
    
    J <- names(rgcca$a)

    title <- ifelse(type == "cor",
            "Variable correlations with",
            "Variable weights on")

    # sort in decreasing order
    df <- data.frame(order_df(df, 1, TRUE), order = NROW(df):1)

    # max threshold for n
    if (NROW(df) >= n_mark)
        df <- df[seq(n_mark), ]

    # if the superblock is selected, color the text of the y-axis according
    # to their belonging to each blocks
    if ((rgcca$superblock && i_block == length(rgcca$a)) || collapse) {
        color <- factor(df$resp)
        levels(color) <- color_group(color)
        p <- ggplot(df, aes(order, df[, 1], fill = df$resp))
    } else {
        color <- "black"
        p <- ggplot(df, aes(order, df[, 1], fill = abs(df[, 1])))
    }

    p <- plot_histogram(
        p,
        df,
        title,
        as.character(color),
        cex = cex,
        cex_sub = cex_sub,
        cex_axis = cex_axis
    ) +
    labs(subtitle = print_comp(rgcca, comp, i_block))

    # If some blocks have any variables in the top hit, selects the ones
    # corresponding
    if (collapse)
        col <- J
    else
        col <- J[-length(J)]

    matched <- match(rev(unique(df$resp)), col)

    # Force all the block names to appear on the legend
    if (length(color) != 1)
        p <- order_color(rgcca$a, p, matched, collapse)

    if ((!rgcca$superblock || i_block != length(rgcca$a)) && !collapse)
            p <- p + theme(legend.position = "none")

    return(p)
}
