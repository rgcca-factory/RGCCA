#' Print bootstrap
#'
#' Print a bootstrap object
#' @param x A fitted bootstrap object (see \code{\link[RGCCA]{bootstrap}})
#' @param type Character string indicating the bootstrapped object to print:
#' block-weight vectors ("weight", default) or block-loading vectors
#' ("loadings").
#' @param empirical A logical value indicating if the bootstrap confidence
#' intervals and p-values are derived from the empirical distribution.
#' (defaut: TRUE)
#' @param display_order A logical value for ordering the variables
#' @param ... Further arguments in print
#' the means, 95\% intervals, bootstrap ratio, p-values and other statistics.
#' @return A matrix containing for each variable:
#' \itemize{
#' \item 'estimate' for block weight or loading vectors.
#' \item 'mean' for the mean of the bootstrap block weight/loading vectors.
#' \item 'sd' for the boostrapped estimate of the standard deviation of the
#' block weight/loading vectors.
#' \item 'lower/upper_bound' for the lower and upper intervals.
#' \item 'pvals' for p-values. In the case of SGCCA, the occurrences of the
#' bootstrap weights are distributed according to the binomial distribution
#' while for RGCCA, their absolute value divided by their standard deviation
#' follows a normal distribution.
#' \item 'adjust.pval' for ajusted p-value (fdr correction by default)
#' }
#' @examples
#' data("Russett")
#' blocks <- list(
#'   agriculture = Russett[, seq(3)],
#'   industry = Russett[, 4:5],
#'   politic = Russett[, 6:11]
#' )
#' fit.rgcca <- rgcca(blocks, ncomp = c(2, 1, 2))
#' boot.out <- bootstrap(fit.rgcca, n_boot = 20, n_cores = 2)
#' print(boot.out)
#' @export
print.bootstrap <- function(x, type = "weight", empirical = TRUE,
                            display_order = FALSE, ...) {
  print_call(x$rgcca$call)
  cat("\n")
  cat(paste0(
    "Extracted statistics on the block-", type, " vectors from ",
    NCOL(x$bootstrap[[1]][[1]][[1]]), " bootstrap samples"
  ), "\n")

  # Remove superblock from the print
  J <- length(x$rgcca$call$raw)
  ncompmax <- max(x$rgcca$call$ncomp[-(J + 1)])

  for (comp in seq(ncompmax)) {
    cat(paste("Component:", comp, "\n"))
    # Extract the blocks for which component comp was extracted
    blocks <- which(vapply(
      x$bootstrap$W[[comp]], function(y) !all(is.na(y)),
      FUN.VALUE = logical(1L)
    ))
    print(Reduce(rbind, lapply(
      blocks,
      function(block) {
        get_bootstrap(
          b = x, type = type,
          block = block,
          comp = comp,
          empirical = empirical,
          display_order = display_order
        )
      }
    )))
  }
}
