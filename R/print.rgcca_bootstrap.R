#' @export
#' @rdname print
#' @order 4
print.rgcca_bootstrap <- function(x, ...) {
  stopifnot(is(x, "rgcca_bootstrap"))
  cat("RGCCA bootstrap model fitted on", x$n_boot, "bootstrap samples.")
  cat("\n")
}
