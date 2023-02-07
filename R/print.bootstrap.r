#' Print rgcca_bootstrap
#'
#' Print a rgcca_bootstrap object
#' @inheritParams plot.bootstrap
#' @param x A fitted rgcca_bootstrap object
#' (see \code{\link[RGCCA]{rgcca_bootstrap}})
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
#' boot.out <- rgcca_bootstrap(fit.rgcca, n_boot = 20, n_cores = 2)
#' print(boot.out)
#' @export
print.bootstrap <- function(x, block = seq_along(x$rgcca$call$blocks),
                            comp = 1, type = "weights", empirical = TRUE,
                            display_order = FALSE, adj.method = "fdr", ...) {
  ### Perform checks and parse arguments
  stopifnot(is(x, "bootstrap"))
  type <- match.arg(type, c("weights", "loadings"))
  lapply(block, function(i) check_blockx("block", i, x$rgcca$call$blocks))
  Map(function(y, z) check_compx(y, y, x$rgcca$call$ncomp, z), comp, block)

  ### Construct data.frame
  column_names <- columns <- c(
    "estimate", "mean", "sd", "lower_bound",
    "upper_bound", "bootstrap_ratio", "pval", "adjust.pval"
  )
  if (!empirical) {
    columns <- c(
      "estimate", "mean", "sd", "th_lower_bound",
      "th_upper_bound", "bootstrap_ratio", "th_pval", "adjust.pval"
    )
  }
  df <- x$stats[x$stats$type == type, ]
  col_pval <- ifelse(empirical, "pval", "th_pval")
  df["adjust.pval"] <- p.adjust(df[, col_pval], method = adj.method)
  df <- df[df$block %in% names(x$rgcca$blocks)[block], ]
  df <- df[df$comp == comp, ]
  rownames(df) <- df$var
  df <- df[, columns]
  colnames(df) <- column_names

  if (display_order) {
    df <- df[order(abs(df$estimate), decreasing = TRUE), ]
  } else {
    df <- df[unlist(lapply(x$rgcca$blocks[block], colnames)), ]
  }

  df <- format(df, digits = 3)

  ### Print
  print_call(x$rgcca$call)
  cat("\n")
  type_str <- ifelse(type == "weights", "weight", "loading")
  cat(paste0(
    "Extracted statistics on the block-", type_str, " vectors for component ",
    comp, " from ", x$n_boot, " bootstrap samples"
  ), "\n")

  print(df, quote = FALSE, ...)
}
