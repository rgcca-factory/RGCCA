#' plot.df_bootstrap
#' @param type "1D" or "2D"
#' @param x result of get_bootstrap
#' @inheritParams plot_bootstrap_1D
#' @inheritParams plot_bootstrap_2D
#' @export
plot.bootstrap=function(x,type="1D",i_block=length(x$rgcca$call$blocks),colors=NULL,title=NULL,cex=1,n_cores= parallel::detectCores() - 1,...)
{
    if(type=="1D")
    {
        p1=plot_bootstrap_1D(
            b = x,
            df_b = NULL,
            x = "estimate",
            y = "occurrences",
            n_mark = 50,
            title = title, 
            colors = colors,
            comp = 1,
            i_block = i_block,
            collapse = FALSE,
            n_cores =n_cores,
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
            cex_main = 25 * cex,
            cex_sub = 16 * cex,
            cex_point = 3 * cex,
            cex_lab = 19 * cex,
            comp = 1,
            i_block = i_block,
            collapse = FALSE,
            n_cores = n_cores)
    }
    plot(p1)
}