#' plot.predict
#' @param x A rgcca_predict object (see \code{\link[RGCCA]{rgcca_predict}})
#' @inheritParams plot_ind
#' @examples
#' data(Russett)
#' @export
plot.predict <- function(x, ...) {
  p <- plot_ind(x$rgcca_res, predicted = x, ...)
  plot(p)
}
