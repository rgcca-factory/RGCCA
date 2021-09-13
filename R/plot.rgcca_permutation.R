#' Plot fitted rgcca permutation object
#'
#' Plot a fitted rgcca permutation object. The various set of tuning parameters
#' are represented in the x-axis and the the RGCCA objective function - obtained
#' from the both the orginal and permuted blocks - in the y-axis. If type =
#' "zstat" the value of the zstat for the various combinations are reported in
#' the y-axis.
#' @param x a fitted rgcca_permutation object
#' (see \code{\link[RGCCA]{rgcca_permutation}})
#' @inheritParams plot_permut_2D
#' @inheritParams plot.rgcca
#' @examples
#' data(Russett)
#' A <- list(agriculture = Russett[, seq(3)],
#'           industry = Russett[, 4:5],
#'           politic = Russett[, 6:11])
#'
#' perm.out <- rgcca_permutation(A, par_type= "tau", n_perms = 2, n_cores = 1)
#' print(perm.out)
#' plot(perm.out)
#'
#' perm.out <- rgcca_permutation(A, par_type = "sparsity",
#'                               n_perms = 5, n_cores = 1)
#' print(perm.out)
#' plot(perm.out, type = "zstat")
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
