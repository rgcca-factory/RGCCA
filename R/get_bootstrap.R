#' Extract statistics from a fitted bootstrap object
#'
#' This function extracts statistical information from a fitted bootstrap
#' object (see \code{\link[RGCCA]{bootstrap}}).
#'
#' @inheritParams bootstrap
#' @inheritParams plot_histogram
#' @inheritParams plot_var_2D
#' @inheritParams plot.rgcca
#' @inheritParams plot_var_1D
#' @param b A fitted bootstrap object (see  \code{\link[RGCCA]{bootstrap}})
#' @param type Character string indicating the bootstrapped object to print:
#' block-weight vectors ("weight", default) or block-loading vectors
#' ("loadings").
#' @param display_order A logical value for ordering the variables
#' @param empirical A logical value indicating if the bootstrap confidence
#' intervals and p-values are derived from the empirical distribution.
#' (defaut: TRUE)
#' @param adj.method character string indicating the method used for p-value
#' adjustment (default: fdr - a.k.a Benjamini-Hochberg correction)
#' @return A dataframe containing:
#' \itemize{
#' \item 'mean' for the mean of the bootstrap weights/loadings
#' (non-null for SGCCA)
#' \item 'estimate' for RGCCA block-weight/block-loading vectors
#' \item 'sd' for the bootstrap estimate of the standard error of the
#' (non-null in case of SGCCA) bootstrap weights/loadings
#' \item 'lower/upper_bound' for the lower and upper intervals
#' (0.025/0.975 percentile).
#' \item 'bootstrap_ratio' defined as the ratio between weight (or loadings)
#' and the bootstrap estimate of the standard deviation.
#' \item 'pval' for p-value (see details)
#' \item 'adjust.pval' for adjusted p-value (default value: fdr (Benjamini-
#' Hochberg correction))
#' \item 'occurrences' for non-zero occurences (for SGCCA)
#' \item 'sign' for significant ('*') or not ('ns') p-value (alpha = 0.05)
#' (see details)
#' }
#' @details
#' For RGCCA, the p-value is computed by assuming that the ratio of the blocks
#' weight values to the bootstrap estimate of the standard deviation follows
#' the standardized normal distribution.
#' By including sparsity (with "sgcca","spls" or "spca"), the frequency of a
#' selected variable may depend on both the level of sparsity and the total
#' number of variables in each block.
#' For a random selection of the variable within the block, the number of
#' occurrences (0 or 1) follows a Bernoulli distribution with the parameter
#' p = proportion of selected variables in the block.
#' This proportion is estimated by the average number of selected variables
#' over all bootstraps divided by the total number of variables in each block
#' (p_j). On a larger number of bootstrap samples, the number of occurrences
#' follows a binomial distribution B(n,p) with n=number of bootstraps.
#' The test is based on the following null hypothesis: "the variable is
#' randomly selected according to B(n,p)". This hypothesis is rejected when
#' the number of occurrences is higher than the 1-(0.05/p_j)th quantile
#' @examples
#' # Bootstrap confidence intervals and p-values for RGCCA
#' data("Russett")
#' blocks <- list(
#'   agriculture = Russett[, seq(3)],
#'   industry = Russett[, 4:5],
#'   politic = Russett[, 6:11]
#' )
#'
#' rgcca_out <- rgcca(blocks)
#' boot <- bootstrap(rgcca_out, n_boot = 5, n_cores = 1)
#' get_bootstrap(boot, type = "loadings")
#'
#' # Stability of the selected variables for SGCCA
#' @export
#' @importFrom stats pt pbinom
#' @seealso \code{\link[RGCCA]{bootstrap}},
#' \code{\link[RGCCA]{plot.bootstrap}},
#' \code{\link[RGCCA]{print.bootstrap}}
get_bootstrap <- function(b, type = "weight", comp = 1,
                          block = length(b$bootstrap[[1]][[1]]),
                          empirical = TRUE,
                          display_order = TRUE,
                          adj.method = "fdr") {
  stopifnot(is(b, "bootstrap"))
  check_ncol(b$rgcca$Y, block)
  check_blockx("block", block, b$rgcca$call$blocks)
  check_compx("comp", comp, b$rgcca$call$ncomp, block)

  if (type == "weight") {
    b$bootstrap <- b$bootstrap$W
  }
  if (type == "loadings") {
    b$bootstrap <- b$bootstrap$L
  }

  bootstrapped <- b$bootstrap[[comp]][[block]]
  n_boot <- sum(!is.na(bootstrapped[1, ]))
  mean <- apply(bootstrapped, 1, function(x) mean(x, na.rm = TRUE))

  if (type == "weight") {
    weight <- b$rgcca$a[[block]][, comp]
  }
  if (type == "loadings") {
    weight <- drop(cor(b$rgcca$call$blocks[[block]],
      b$rgcca$Y[[block]][, comp],
      use = "pairwise.complete.obs"
    ))
  }

  if (empirical) {
    if (type == "weight") {
      std <- apply(bootstrapped, 1, function(x) sd(x, na.rm = TRUE))
      bootstrap_ratio <- weight / std
      cross_zero <- apply(
        bootstrapped, 1,
        function(x) {
          c(
            sum(x > 0, na.rm = TRUE),
            sum(x < 0, na.rm = TRUE)
          )
        }
      )
      p.vals <- apply(cross_zero, 2, function(x) min(x) / max(x))
      lower_bound <- apply(
        bootstrapped, 1,
        function(y) {
          return(quantile(y, 0.025))
        }
      )

      upper_bound <- apply(
        bootstrapped, 1,
        function(y) {
          return(quantile(y, 0.975))
        }
      )
    }

    if (type == "loadings") {
      ftrans <- 0.5 * log((1 + weight) / (1 - weight))
      r <- 0.5 * log((1 + bootstrapped) / (1 - bootstrapped))
      std <- apply(r, 1, function(x) sd(x, na.rm = TRUE))
      bootstrap_ratio <- ftrans / std
      cross_zero <- apply(
        bootstrapped, 1,
        function(x) {
          c(
            sum(x > 0, na.rm = TRUE),
            sum(x < 0, na.rm = TRUE)
          )
        }
      )
      p.vals <- apply(cross_zero, 2, function(x) min(x) / max(x))

      lower_bound <- apply(
        bootstrapped, 1,
        function(y) {
          return(quantile(y, 0.025))
        }
      )

      upper_bound <- apply(
        bootstrapped, 1,
        function(y) {
          return(quantile(y, 0.975))
        }
      )
    }
  } else {
    tail <- qnorm(1 - .05 / 2)

    if (type == "weight") {
      std <- apply(bootstrapped, 1, function(x) sd(x, na.rm = TRUE))
      bootstrap_ratio <- weight / std
      p.vals <- 2 * pnorm(abs(bootstrap_ratio), lower.tail = FALSE)
      lower_bound <- weight - std * tail
      upper_bound <- weight + std * tail
    }

    if (type == "loadings") {
      ftrans <- 0.5 * log((1 + weight) / (1 - weight))
      r <- 0.5 * log((1 + bootstrapped) / (1 - bootstrapped))
      std <- apply(r, 1, function(x) sd(x, na.rm = TRUE))
      bootstrap_ratio <- ftrans / std
      p.vals <- 2 * pnorm(abs(bootstrap_ratio), lower.tail = FALSE)
      lower_bound <- weight - std * tail
      upper_bound <- weight + std * tail
    }
  }

  df <- data.frame(
    estimate = weight,
    mean = mean,
    sd = std,
    lower_bound = lower_bound,
    upper_bound = upper_bound,
    bootstrap_ratio = bootstrap_ratio,
    pval = p.vals,
    adjust.pval = p.adjust(p.vals, method = adj.method)
  )

  if (display_order) {
    index <- which(colnames(df) == "estimate")
    df <- data.frame(order_df(df, index, allCol = TRUE))
  }

  attributes(df)$indexes <- list(
    estimate = ifelse(type == "weight",
      "block-weight",
      "block-loadings"
    ),
    bootstrap_ratio = "Bootstrap-ratio",
    sign = "Significance",
    mean = "Mean bootstrap weights"
  )

  attributes(df)$method <- class(b$rgcca)
  attributes(df)$n_boot <- n_boot
  attributes(df)$n_blocks <- length(b$rgcca$call$blocks)
  class(df) <- c(class(df), "df_bootstrap")
  return(df)
}
