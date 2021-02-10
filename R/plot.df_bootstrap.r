#' Plot a fitted bootstrap object
#' 
#' Plot the results of a fitted bootstrap object. 
#' @inheritParams plot_var_2D
#' @inheritParams plot_var_1D
#' @inheritParams plot.rgcca
#' @param x A fitted bootstrap object (see \code{\link[RGCCA]{bootstrap}})
#' @inheritParams plot_bootstrap_1D
#' @inheritParams plot_bootstrap_2D
#' @inheritParams plot.rgcca
#' @export
#' @examples
#' data("Russett")
#' blocks = list(agriculture = Russett[, seq(3)], 
#'               industry = Russett[, 4:5], 
#'               politic = Russett[, 6:11])
#' fit.rgcca = rgcca(blocks)
#' fit.boot=bootstrap(fit.rgcca, n_boot = 2, n_cores = 1)
#' plot(fit.boot)

plot.bootstrap=function(x, block = length(x$rgcca$call$blocks), 
                        comp = 1, n_mark = 30, bars = "quantile",
                        colors = NULL, title = NULL, cex = 1,
                        collapse = FALSE, cex_main = 14, cex_sub = 12, 
                        cex_point = 10, cex_lab = 10, cex_axis = 10, ...)
{
    stopifnot(is(x, "bootstrap"))
    check_blockx("block", block, x$rgcca$call$blocks)
    
    if(x$rgcca$call$superblock)
    {
        if(block==length(x$rgcca$call$blocks))
            block=length(x$rgcca$call$blocks)-1
    }

    if(x$rgcca$call$type%in%c("sgcca", "spls", "spca"))
        {
            x1="occurrences"
            y1="estimate"
            title=ifelse(is.null(title),
                         paste0("Non-zero occurrences (", 
                                names(x$rgcca$call$blocks)[block], 
                                ")\n(", ncol(x$bootstrap[[1]][[1]]),
                                " bootstrap samples - comp ", comp, ")"), 
                         title)
        }
         else
         {
            x1="estimate"
            y1="sign"
            title=ifelse(is.null(title), 
                         paste0("Bootstrap confidence interval (",
                                names(x$rgcca$call$blocks)[block], ")\n (", 
                                ncol(x$bootstrap[[1]][[1]]),
                                " bootstrap samples - comp ", comp, ")" ), 
                         title)
         }
    
           
         p1 = plot_bootstrap_1D(
              b = x,
              df_b = NULL,
              x = x1,
              y = y1,
              n_mark = n_mark,
              title = title, 
              colors = colors,
              comp = comp,
              i_block = block,
              collapse = collapse,
              bars = bars,
              cex = 1,
              cex_main = cex_main,
              cex_sub = cex_sub,
              cex_axis = cex_axis)
  
    plot(p1)
}
