#' Plot permutation
#' 
#' Plot a permutation object (tuning RGCCA parameters). The parameters tuned 
#' for maximizing RGCCA criteria is displayed in the title. In abscissa, the
#' index of combination (number corresponding to the rownames of 'penalties' 
#' in the tuning parameters object). In ordinate, a score depending of the
#' 'type' parameter. For 'type = crit', the grey dots correspond to the RGCCA 
#' criteria of each permutation for each combination. For 'type = zstats', 
#' the three horizontal bars in dashed lines correspond to the three 
#' significance thresholds from the lowest  to the highest bar: 95%, 99% and 99.9%.
#'  The best parameters are in red by default.
#' @param x result of rgcca_permutation (see  \code{\link[RGCCA]{rgcca_permutation}})
#' @inheritParams plot_permut_2D
#' @inheritParams plot.rgcca
#' @examples
#' data("Russett")
#' A <- list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
#'     politic = Russett[, 6:11] )
#' perm <- rgcca_permutation(A,par_type="tau", n_run = 2, n_cores = 1)
#' plot(perm)
#' perm <- rgcca_permutation(A, par_type = "sparsity", n_run = 5, n_cores = 1)
#' plot(perm, type="zstat")
#' plot(perm)
#' @export
#' @importFrom ggplot2 ggplot
plot.permutation = function(x,
    type = "crit",
    cex = 1,
    title = NULL,
    cex_main = 14 * cex,
    cex_sub = 12 * cex,
    bars = "points" ,
    cex_point = 3 * cex,
    cex_lab = 19 * cex,
    colors = c("red", "grey"),
    ...){

    p1 <- plot_permut_2D(
        x,
        type = type,
        cex = cex,
        cex_main = cex_main,
        cex_sub = cex_sub,
        bars = bars,
        cex_point = cex_point,
        cex_lab = cex_lab,
        colors  = colors,
        title = title
    )
    plot(p1)
}