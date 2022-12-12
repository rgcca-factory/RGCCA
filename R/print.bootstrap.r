#' Print bootstrap
#'
#' Print a bootstrap object
#' @param x A fitted bootstrap object (see \code{\link[RGCCA]{bootstrap}})
#' @param type Character string indicating the bootstrapped object to print:
#' block-weight vectors ("weights", default) or block-loading vectors
#' ("loadings").
#' @param empirical A logical value indicating if the bootstrap confidence
#' intervals and p-values are derived from the empirical distribution.
#' (defaut: TRUE)
#' @param display_order A logical value for ordering the variables
#' @param ... Further arguments in print
#' the means, 95\% intervals, bootstrap ratio, p-values and other statistics.
#' @return none
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
print.bootstrap <- function(x, type = "weights", empirical = TRUE,
                            display_order = FALSE, ...) {
  type <- match.arg(type, c("weights", "loadings"))
  print_call(x$rgcca$call)
  cat("\n")
  type_str <- ifelse(type == "weights", "weight", "loading")
  cat(paste0(
    "Extracted statistics on the block-", type_str, " vectors from ",
    NCOL(x$bootstrap[[1]][[1]][[1]]), " bootstrap samples"
  ), "\n")

  # Remove superblock from the print
  J <- length(x$rgcca$call$blocks)
  ncompmax <- max(x$rgcca$call$ncomp[-(J + 1)])

  for (comp in seq(ncompmax)) {
    cat(paste("Component:", comp, "\n"))
    # Extract the blocks for which component comp was extracted
    blocks <- which(vapply(
      x$bootstrap$W[[comp]], function(y) !all(is.na(y)),
      FUN.VALUE = logical(1L)
    ))

    df <- Reduce(rbind, lapply(
      blocks,
      function(block) {
        get_bootstrap(
          b = x, type = type,
          block = block,
          comp = comp,
          empirical = empirical
        )
      }
    ))
    if (display_order) {
      df <- df[order(abs(df$estimate), decreasing = TRUE), ]
    }
    print(df)
  }
}
