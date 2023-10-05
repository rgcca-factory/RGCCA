#' @export
#' @rdname summary
#' @order 5
summary.rgcca_stability <- function(object, ...) {
  stopifnot(is(object, "rgcca_stability"))
  print(object$rgcca_res, ...)
  cat("\n")
}
