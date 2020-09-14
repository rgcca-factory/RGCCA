#' plot.cv
#' 
#' Plots the results of rgcca_cv (in abscissa, the configurations, in ordinates the RMSE on test dataset)
#' @param x result of rgcca_crossvalidation (see  \code{\link[RGCCA]{rgcca_crossvalidation}} )
#' @inheritParams plot_ind
#' @examples
#' data(Russett)
#'blocks = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
#'              politic = Russett[, 6:11] )
#' rgcca_out = rgcca(blocks, response = 3,ncomp=2)
#' cv=rgcca_crossvalidation(rgcca_out, validation = "kfold", k = 5, n_cores = 1)
#' plot(cv)
#' @export
plot.cv=function(x, ...) {
       plot_ind(x$rgcca_res,predicted=x,...)
  # if(type=="error")
  # {
  #     l=list(rgcca0=x$list_rgcca[[1]],rgccaList=x$list_rgcca )
  #     class(l)="list_rgcca"
  #     plot(l,type="ind",compx=1,compy=1)
  # }
}