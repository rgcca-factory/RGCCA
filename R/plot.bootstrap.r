#' Plot a fitted bootstrap object
#'
#' Plot the results of a fitted bootstrap object.
#' @inheritParams plot_var_2D
#' @inheritParams plot_var_1D
#' @inheritParams plot.rgcca
#' @inheritParams plot2D
#' @param x A fitted bootstrap object (see \code{\link[RGCCA]{bootstrap}})
#' @param type Character string indicating the bootstrapped object to plot: 
#' block-weight vectors ("weight", default) or block-loading vectors 
#' ("loadings"). 
#' @export
#' @examples
#' data("Russett")
#' blocks = list(agriculture = Russett[, seq(3)],
#'               industry = Russett[, 4:5],
#'               politic = Russett[, 6:11])
#' fit.rgcca = rgcca(blocks, ncomp = 2, method = "rgcca", tau = 1)
#' fit.boot=bootstrap(fit.rgcca, n_boot = 20, n_cores = 1)
#' plot(fit.boot, type = "weight", block = 1, comp=1)

plot.bootstrap=function(x, block = length(x$rgcca$call$blocks),
                        comp = 1, type = "weight", n_mark = 30, 
                        colors = NULL, title = NULL, cex = 1, 
                        cex_main = 14, cex_sub = 12,
                        cex_point = 10, cex_lab = 10, cex_axis = 10, ...)
{
    stopifnot(is(x, "bootstrap"))
    check_blockx("block", block, x$rgcca$call$blocks)

    if(x$rgcca$call$superblock)
    {
        if(block==length(x$rgcca$call$blocks))
            block=length(x$rgcca$call$blocks)-1
    }

    x1="estimate"
    y1="sign"
    title=ifelse(is.null(title),
                 paste0("Bootstrap confidence interval (",
                         names(x$rgcca$call$blocks)[block], ")\n (",
                         type, " - ",
                         ncol(x$bootstrap[[1]][[1]][[1]]),
                         " bootstrap samples - comp ", comp, ")" ),
                         title)

     p1 = plot_bootstrap_1D(
              b = x,
              df_b = NULL,
              type = type, 
              x = x1,
              y = y1,
              n_mark = n_mark,
              title = title,
              colors = colors,
              comp = comp,
              i_block = block, 
              cex = 1,
              cex_main = cex_main,
              cex_sub = cex_sub,
              cex_axis = cex_axis)

    plot(p1)
}
