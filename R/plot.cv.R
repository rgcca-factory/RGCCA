#' plot.cv
#' @param x result of rgcca_crossvalidation (see  \code{\link[RGCCA]{rgcca_crossvalidation}} )
#' @inheritParams plot_ind
#' @examples
#' data(Russett)
#'blocks = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
#'              politic = Russett[, 6:11] )
#' rgcca_out = rgcca(blocks, response = 3)
#' cv=rgcca_crossvalidation(rgcca_out, validation = "kfold", k = 5, n_cores = 1)
#' plot(cv)
#' @export
plot.cv=function(x,...)
{
   
   p<- plot_ind(x$rgcca_res,predicted=x,...)
   plot(p)
 
}