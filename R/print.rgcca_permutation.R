#' @export
#' @rdname print
#' @order 3
print.rgcca_permutation <- function(x, ...) {
  stopifnot(is(x, "rgcca_permutation"))
  cat(
    "RGCCA permutation model fitted on", nrow(x$params),
    "sets of parameters with", x$n_perms, "permutations each."
  )
  cat("\n")
}
