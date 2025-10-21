#' Compute bootstrap stats
#'
#' This function aims to extract statistics over the different bootstrap
#' samples.
#' @param res A data.frame of raw bootstrap results.
#' @param rgcca_res A fitted rgcca object.
#' @param n_boot Number of bootstrap samples.
#' @return A data.frame containing the computed statistics.
#' @importFrom matrixStats rowSums2 rowMeans2 rowSds
#' @importFrom matrixStats rowMins rowMaxs rowQuantiles
#' @noRd
rgcca_bootstrap_stats <- function(res, rgcca_res, n_boot) {
  ### Compute aggregated statistics
  # Utility function to compute the statistics
  aggregated_stats <- function(df, type) {
    if (type == "weights") {
      y <- df$value
    } else {
      y <- 0.5 * log((1 + df$value) / (1 - df$value))
    }
    count_pos <- rowSums2(matrix(df$value > 0, ncol = n_boot))
    count_neg <- rowSums2(matrix(df$value < 0, ncol = n_boot))
    z <- cbind(count_pos, count_neg)
    res <- cbind(
      rowMeans2(matrix(df$value, ncol = n_boot)),
      rowSds(matrix(y, ncol = n_boot)),
      rowMins(z) / rowMaxs(z),
      rowQuantiles(
        matrix(df$value, ncol = n_boot),
        probs = c(0.025, 0.975)
      )
    )
    colnames(res) <- c("mean", "sd", "pval", "lower_bound", "upper_bound")

    res <- data.frame(res)
    res$var <- matrix(df$var, ncol = n_boot)[, 1]
    res$comp <- matrix(df$comp, ncol = n_boot)[, 1]
    res$block <- matrix(df$block, ncol = n_boot)[, 1]

    res$type <- type

    return(res)
  }

  stats <- rbind(
    aggregated_stats(res[res$type == "weights", ], "weights"),
    aggregated_stats(res[res$type == "loadings", ], "loadings")
  )

  ### Compute additional statistics using aggregated statistics
  tail <- qnorm(1 - .05 / 2)
  stats$estimate <- c(
    unlist(rgcca_res$a, use.names = FALSE),
    unlist(lapply(unique(stats$block), function(block) {
      cor2(
        to_mat(rgcca_res$blocks[[block]]),
        rgcca_res$Y[[block]]
      )
    }), use.names = FALSE)
  )
  suppressWarnings(
    ftrans <- 0.5 * log((1 + stats$estimate) / (1 - stats$estimate))
  )
  ftrans[is.na(ftrans)] <- 0
  ftrans <- ftrans * (stats$type == "loadings") +
    stats$estimate * (stats$type == "weights")
  stats$bootstrap_ratio <- ftrans / stats$sd
  stats$th_pval <- 2 * pnorm(abs(stats$bootstrap_ratio), lower.tail = FALSE)
  stats$th_lower_bound <- stats$estimate - stats$sd * tail
  stats$th_upper_bound <- stats$estimate + stats$sd * tail

  return(stats)
}
