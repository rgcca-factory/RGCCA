#' Barplot of a fingerprint
#'
#' Barplot of the higher outer weight vectors for a component of a block 
#' (by default, the superblock or the last one) analysed by R/SGCCA
#'
#' @inheritParams plot_var_2D
#' @inheritParams plot_histogram
#' @param comp An integer giving the index of the analysis components
#' @param type A string giving the criterion to selects biomarkers : either 
#' "cor" for correlation between the component and the block
#' or "weight" for the weight of the RGCCA
#' @seealso \code{\link[RGCCA]{rgccad}}, \code{\link[RGCCA]{sgcca}}
#' @examples
#' weights = lapply(seq(3), function(x) matrix(runif(7*2), 7, 2))
#' for (i in seq(3))
#' row.names(weights[[i]]) <- paste0(letters[i],
#'      letters[seq(NROW(weights[[i]]))])
#' weights[[4]] = Reduce(rbind, weights)
#' rgcca_out = list(a = weights, call = list(type="rgcca", ncomp = rep(2,4)))
#' names(rgcca_out$a) = LETTERS[seq(4)]
#' rgcca_out$call$blocks = lapply(rgcca_out$a, t)
#' rgcca_out$call$superblock = TRUE
#' class(rgcca_out) = "rgcca"
#' # With the 1rst component of the superblock
#' plot_var_1D(rgcca_out, 1, type = "weight")
#' # With the 2nd component of the 1rst block by selecting the ten higher weights
#' plot_var_1D(rgcca_out, 2, 10, 1, type = "weight")
#' data("Russett")
#' blocks = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
#'     politic = Russett[, 6:11] )
#' rgcca_out = rgcca(blocks)
#' plot_var_1D(rgcca_out, collapse = TRUE)
#' @export
#' @importFrom ggplot2 ggplot
plot_var_1D <- function(
    rgcca_res,
    comp = 1,
    n_mark = 100,
    i_block = length(rgcca_res$a),
    type = "cor",
    collapse = FALSE,
    title = NULL,
    colors = NULL,
    ...) {

    check_colors(colors)
    if(rgcca_res$call$superblock==FALSE){collapse=FALSE}
    
    df <- get_ctr2(
        rgcca_res = rgcca_res,
        compx = comp,
        compy = comp,
        i_block = i_block,
        type = type,
        n_mark = n_mark,
        collapse = collapse,
        remove_var = FALSE
    )
    resp <- df$resp

    if (i_block < length(rgcca_res$a) || tolower(rgcca_res$call$type) == "pca")
        rgcca_res$call$superblock <- FALSE
     J <- names(rgcca_res$a)

    if (is.null(title))
        title <- ifelse(type == "cor",
            "Variable correlations",
            "Variable weights")

    # sort in decreasing order
    df <- data.frame(order_df(df, 1, TRUE), order = NROW(df):1)
    class(df) <- c(class(df), "d_var1D")

    # max threshold for n
    if (NROW(df) >= n_mark)
        df <- df[seq(n_mark), ]

    # if the superblock is selected, color the text of the y-axis according
    # to their belonging to each blocks

    if ((rgcca_res$call$superblock && i_block == length(rgcca_res$a)) || collapse) {
        color <- factor(resp)
        levels(color) <- color_group(color, colors = colors)
        p <- ggplot(df, aes(order, df[, 1], fill = resp))
    } else {
        color <- "black"
        p <- ggplot(df, aes(order, df[, 1], fill = abs(df[, 1])))
    }

    p <- plot_histogram(
        p,
        df,
        title,
        as.character(color),
        colors = colors,
        ...
    ) +
    labs(subtitle = print_comp(rgcca_res, comp, i_block))

    # If some blocks have any variables in the top hit, selects the ones
    # corresponding
    if (collapse)
        col <- J
    else
        col <- J[-length(J)]

    matched <- match(rev(unique(resp)), col)

    # Force all the block names to appear on the legend
    if (length(color) != 1)
        p <- suppressMessages(order_color(rgcca_res$a, p, matched, collapse, colors))

    if ((!rgcca_res$call$superblock || i_block != length(rgcca_res$a)) && !collapse)
            p <- p + theme(legend.position = "none")

    return(p)
}
