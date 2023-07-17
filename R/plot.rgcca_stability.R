#' @export
#' @rdname plot
#' @order 5
plot.rgcca_stability <- function(x, ...) {
  stopifnot(is(x, "rgcca_stability"))
  plot(x$rgcca_res, ...)
}
