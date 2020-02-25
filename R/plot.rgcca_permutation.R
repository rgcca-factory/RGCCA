#' plot.permutation
#' @param x result of rgcca_permutation
#' @param main title of the plot
#' @inheritParams plot_permut_2D
#' @param ... Further graphical parameters
#' @examples
#' data("Russett")
#' A = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
#'     politic = Russett[, 6:11] )
#' perm <- rgcca_permutation(A, nperm = 2, n_cores = 1)
#' plot(perm)
#' perm <- rgcca_permutation(A, p_spars = TRUE, nperm = 2, n_cores = 1)
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
            ...
        )
        plot(p1)
        invisible(p1)
   
}