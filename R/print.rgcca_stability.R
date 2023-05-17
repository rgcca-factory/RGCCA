#' @export
#' @rdname print
#' @order 5
print.rgcca_stability <- function(x, ...) {
  stopifnot(is(x, "rgcca_stability"))
  print(x$rgcca_res, ...)
}
