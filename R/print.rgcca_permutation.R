#' @export
#' @rdname print
#' @order 3
print.rgcca_permutation <- function(x, ...) {
  stopifnot(is(x, "rgcca_permutation"))
  cat(
    "RGCCA permutation object obtained with", nrow(x$params),
    "sets of parameters and", x$n_perms, "permutations each."
  )
  cat("\n")
}
