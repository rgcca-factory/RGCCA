#' Plot a bootstrap object
#' 
#' Plot the results of a bootstrap object. The representation can be with bars
#' (1D) or a biplot (2D) (see details).
#' @inheritParams plot_var_2D
#' @inheritParams plot_var_1D
#' @param type A character to select a type of plot among "1D" or "2D"
#' @param x A bootstrap object (see \code{\link[RGCCA]{bootstrap}})
#' @inheritParams plot_bootstrap_1D
#' @inheritParams plot_bootstrap_2D
#' @inheritParams plot.rgcca
#' @export
#' @examples
#' data("Russett")
#' blocks = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
#'     politic = Russett[, 6:11] )
#' rgcca_out = rgcca(blocks)
#' b=bootstrap(rgcca_out, n_boot = 2, n_cores = 1)
#' plot(b,n_cores=1)
#' @return  A visualization of bootstrap results
#' \itemize{
#' \item If type = '1D': barplot of the best variables from a bootstrap with, on the x-axis,
#' the number of non-zero occurrences (SGCCA) or the mean of the bootstrap weights 
#' (RGCCA). The bars are colored according to the significant 95% bootstrap 
#' intervals ('*' or 'ns'; see 'p.vals' in details). In SGCA, the significant variables 
#'  are those above the three bars, respectively, with an alpha = 0.05 
#'  (dark red), 0.01 (red) and 0.001 (light red).
#' \item type = '2D' : graph of the best variables from a bootstrap with, in x-axis, the number of
#' non-zero occurences (SGCCA) or the significant 95% bootstrap 
#' intervals (RGCCA). In in y-axis are the bootstrap-ratios (mean/sd) 
#' Negative weights are colored in red and the positive ones are in green.
#' }
#' @details
#' The function \code{\link[RGCCA]{get_bootstrap}} allows obtaining the numeric values used to produce the graph

plot.bootstrap=function(x,type="1D",block=length(x$rgcca$call$blocks),comp=1,n_mark=30,bars="quantile",colors=NULL,title=NULL,cex=1,n_cores= parallel::detectCores() - 1,...)
{
    stopifnot(is(x, "bootstrap"))
    check_blockx("block", block, x$rgcca$call$blocks)
    match.arg(type,c("1D","2D"))
    if(type=="1D")
    {
        
        if(x$rgcca$call$type%in%c("sgcca","spls","spca"))
        {
            x1="occurrences";
            y1="estimate";
            title=ifelse(is.null(title),
                         paste0("Occurrences: ",names(x$rgcca$call$blocks)[block], " \n(", ncol(x$bootstrap[[1]][[1]])," bootstraps)"),title)
        }
         else
         {
            x1="estimate";
            y1="sign";
            title=ifelse(is.null(title),paste0("Weights: ",names(x$rgcca$call$blocks)[block],"\n (", ncol(x$bootstrap[[1]][[1]])," bootstraps)"),title)}
           
         p1=plot_bootstrap_1D(
            b = x,
            df_b = NULL,
            x = x1,
            y = y1,
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
        if(tolower(x$rgcca$call$type) %in% c("spls", "spca", "sgcca"))
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
         #   sub="Green line indicates significance at 0.05 with Bonferroni correction")
        }
      else
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
        #   sub="Green line indicates significance at 0.05 with Bonferroni correction")
    }
    }
    plot(p1)
}
