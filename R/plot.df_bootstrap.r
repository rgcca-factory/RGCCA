#' Plots a bootstrap object
#' 
#' Plots the results of a bootstrap object. The representation can be with bars (1D) or (2D) #TODO
#' @param type "1D" or "2D"
#' @param x result of bootstrap  \code{\link[RGCCA]{bootstrap}} 
#' @param block number of the block to be plotted
#' @param comp number of the component to be plotted
#' @inheritParams plot_bootstrap_1D
#' @inheritParams plot_bootstrap_2D
#' @export
#' @examples
#' library(RGCCA)
#' data("Russett")
#' blocks = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
#'     politic = Russett[, 6:11] )
#' rgcca_out = rgcca(blocks)
#' b=bootstrap(rgcca_out, n_boot = 2, n_cores = 1)
#' plot(b,n_cores=1)
plot.bootstrap=function(x,type="1D",block=length(x$rgcca$call$blocks),comp=1,n_mark=30,bars="sd",colors=NULL,title=NULL,cex=1,n_cores= parallel::detectCores() - 1,...)
{
    
    if(type=="1D")
    {
            p1=plot_bootstrap_1D(
            b = x,
            df_b = NULL,
            x = "estimate",
            y = "occurrences",
            n_mark = n_mark,
            title = title, 
            colors = colors,
            comp = comp,
            i_block = block,
            collapse = FALSE,
            n_cores =n_cores,
            bars=bars,
            ...)
    }
    if(type=="2D")
    {
        p1=plot_bootstrap_2D(
            b =x,
            df_b = NULL,
            x = "bootstrap_ratio",
            y = "occurrences",
            title = paste("Variable selection \nby",
                          attributes(x)$n_boot,
                          "bootstraps"),
            colors = NULL,
            cex = cex,
            cex_main = 14 * cex,
            cex_sub = 12 * cex,
            cex_point = 10 * cex,
            cex_lab = 10 * cex,
            comp = comp,
            i_block = block,
            collapse = FALSE,
            n_cores = n_cores)
    }
    plot(p1)
}