#' Print a rgcca_stability object.
#'
#' The fitted RGCCA model returned by rgcca_stability is printed
#' All arguments are forwarded to the print.rgcca function.
#' @param x Object of type "stability" produced by rgcca_stability.
#' @param ... Arguments for the print.rgcca function.
#' @return none
#' @examples
#' data(Russett)
#' blocks <- list(
#'   agriculture = Russett[, seq(3)],
#'   industry = Russett[, 4:5],
#'   politic = Russett[, 6:11]
#' )
#' fit.sgcca <- rgcca(blocks, sparsity = c(.8, .9, .6))
#' res <- rgcca_stability(fit.sgcca, n_boot = 10, verbose = FALSE)
#' print(res)
#' @export
print.stability <- function(x, ...) {
  stopifnot(is(x, "stability"))
  print(x$rgcca_res, ...)
}
