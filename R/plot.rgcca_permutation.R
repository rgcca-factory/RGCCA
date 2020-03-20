#' plot.permutation
#' Plots a permutation object. The parameters tuned for maximizing RGCCA criteria is displayed in the title. 
#' In x, the index of combination (number corresponding to the rownames of tuning parameters object). In ordinate, a score depending of the type parameter (RGCCA criterion for crit, and zstat for the pseudo z-scores)
#' @param x result of rgcca_permutation (see  \code{\link[RGCCA]{rgcca_permutation}} )
#' @param main title of the plot
#' @inheritParams plot_permut_2D
#' @param ... Further graphical parameters
#' @examples
#' data("Russett")
#' A = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
#'     politic = Russett[, 6:11] )
#' perm <- rgcca_permutation(A, nperm = 2, n_cores = 1)
#' plot(perm)
#' perm <- rgcca_permutation(A, perm.par = "sparsity", nperm = 5, n_cores = 1)
#' plot(perm,type="crit")
#' plot(perm)
#' @export
plot.permutation=function(x,type="zstat",cex = 1, main= NULL, cex_main = 25 * cex,                          cex_sub = 16 * cex,
                          cex_point = 3 * cex, cex_lab = 19 * cex,...)
{
        p1 <- plot_permut_2D(
            x, 
            type = type,
            cex = 1,
            title = main,
            bars="points",
            ...
        )
        plot(p1)
        invisible(p1)
   
}