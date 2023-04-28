#' Plot a rgcca_stability object.
#'
#' The fitted RGCCA model returned by rgcca_stability is plotted.
#' All arguments are forwarded to the plot.rgcca function.
#' @param x Object of type "stability" produced by rgcca_stability.
#' @param ... Arguments for the plot.rgcca function.
#' @return A ggplot2 plot object.
#' @examples
#' data(Russett)
#' blocks <- list(
#'   agriculture = Russett[, seq(3)],
#'   industry = Russett[, 4:5],
#'   politic = Russett[, 6:11]
#' )
#' fit.sgcca <- rgcca(blocks, sparsity = c(.8, .9, .6))
#' res <- rgcca_stability(
#'   fit.sgcca, n_boot = 10, verbose = TRUE, keep = rep(.1, 3)
#' )
#' plot(res, type = "weights")
#' @export
plot.stability <- function(x, ...) {
  stopifnot(is(x, "stability"))
  plot(x$rgcca_res, ...)
}
