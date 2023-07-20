#' @export
#' @rdname print
#' @order 4
print.rgcca_bootstrap <- function(x, ...) {
  stopifnot(is(x, "rgcca_bootstrap"))
  cat("RGCCA bootstrap object obtained with", x$n_boot, "bootstrap samples.")
  cat("\n")
}
