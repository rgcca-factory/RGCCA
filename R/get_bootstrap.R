#' Extract statistics from a fitted bootstrap object
#'
#' This function extracts statistical information from a fitted bootstrap
#' object (see \code{\link[RGCCA]{bootstrap}}).
#'
#' @inheritParams plot.rgcca
#' @param b A fitted bootstrap object (see  \code{\link[RGCCA]{bootstrap}})
#' @param type Character string indicating the bootstrapped object to print:
#' block-weight vectors ("weights", default) or block-loading vectors
#' ("loadings").
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
#' @seealso \code{\link[RGCCA]{bootstrap}},
#' \code{\link[RGCCA]{plot.bootstrap}},
#' \code{\link[RGCCA]{print.bootstrap}}
get_bootstrap <- function(b, type = "weights", comp = 1,
                          block = 1,
                          empirical = TRUE,
                          adj.method = "fdr") {
  ### Auxiliary functions to compute statistics
  empirical_statistics <- function(x, y, std) {
    bootstrap_ratio <- y / std
    cross_zero <- apply(x, 1, function(z) {
      c(
        sum(z > 0, na.rm = TRUE),
        sum(z < 0, na.rm = TRUE)
      )
    })
    p.vals <- apply(cross_zero, 2, function(z) min(z) / max(z))
    lower_bound <- apply(x, 1, quantile, 0.025)
    upper_bound <- apply(x, 1, quantile, 0.975)

    return(list(
      bootstrap_ratio = bootstrap_ratio,
      p.vals = p.vals,
      lower_bound = lower_bound,
      upper_bound = upper_bound
    ))
  }

  theoretical_statistics <- function(y, std, tail, z = y) {
    bootstrap_ratio <- z / std
    p.vals <- 2 * pnorm(abs(bootstrap_ratio), lower.tail = FALSE)
    lower_bound <- y - std * tail
    upper_bound <- y + std * tail

    return(list(
      bootstrap_ratio = bootstrap_ratio,
      p.vals = p.vals,
      lower_bound = lower_bound,
      upper_bound = upper_bound
    ))
  }

  ### Perform checks
  stopifnot(is(b, "bootstrap"))
  check_blockx("block", block, b$rgcca$call$blocks)
  check_compx("comp", comp, b$rgcca$call$ncomp, block)

  ### Get bootstrap object and estimate
  if (type == "weights") {
    b$bootstrap <- b$bootstrap$W
  } else {
    b$bootstrap <- b$bootstrap$L
  }

  bootstrapped <- b$bootstrap[[comp]][[block]]

  if (type == "weights") {
    estimate <- b$rgcca$a[[block]][, comp]
  } else {
    estimate <- drop(cor(b$rgcca$blocks[[block]],
      b$rgcca$Y[[block]][, comp],
      use = "pairwise.complete.obs"
    ))
  }

  ### Compute statistics
  mean <- apply(bootstrapped, 1, function(x) mean(x, na.rm = TRUE))
  tail <- qnorm(1 - .05 / 2)

  if (type == "weights") {
    std <- apply(bootstrapped, 1, function(x) sd(x, na.rm = TRUE))
    if (empirical) {
      statistics <- empirical_statistics(bootstrapped, estimate, std)
    } else {
      statistics <- theoretical_statistics(estimate, std, tail)
    }
  } else {
    ftrans <- 0.5 * log((1 + estimate) / (1 - estimate))
    r <- 0.5 * log((1 + bootstrapped) / (1 - bootstrapped))
    std <- apply(r, 1, function(x) sd(x, na.rm = TRUE))
    if (empirical) {
      statistics <- empirical_statistics(bootstrapped, ftrans, std)
    } else {
      statistics <- theoretical_statistics(estimate, std, tail, ftrans)
    }
  }

  ### Construct summary data frame
  df <- data.frame(
    estimate = estimate,
    mean = mean,
    sd = std,
    lower_bound = statistics$lower_bound,
    upper_bound = statistics$upper_bound,
    bootstrap_ratio = statistics$bootstrap_ratio,
    pval = statistics$p.vals,
    adjust.pval = p.adjust(statistics$p.vals, method = adj.method)
  )

  return(df)
}
