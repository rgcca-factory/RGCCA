#' Compute bootstrap stats
#'
#' This function aims to extract statistics over the different bootstrap
#' samples.
#' @param res A data.frame of raw bootstrap results.
#' @param rgcca_res A fiited rgcca object.
#' @param verbose A logical indicating if progress bars should be shown.
#' @return A data.frame containing the computed statistics.
#' @noRd
rgcca_bootstrap_stats <- function(res, rgcca_res, verbose) {
  ### Compute aggregated statistics
  # Utility function to compute the statistics
  aggregated_stats <- function(x, type) {
    if (type == "weights") {
      sd <- sd(x)
    } else {
      sd <- sd(0.5 * log((1 + x) / (1 - x)))
    }
    z <- c(
      sum(x > 0, na.rm = TRUE),
      sum(x < 0, na.rm = TRUE)
    )
    list(
      mean = mean(x),
      sd = sd,
      pval = min(z) / max(z),
      lower_bound = quantile(x, 0.025),
      upper_bound = quantile(x, 0.975)
    )
  }

  # Set progress bar option
  if (!verbose) {
    pbapply::pboptions(type = "none")
  } else {
    pbapply::pboptions(type = "timer")
  }

  # Compute var2block to later retrieve "block" from "var"
  var2block <- subset(
    res, res$type == "weights" & res$comp == 1 & res$boot == 1
  )[, c("var", "block")]
  rownames(var2block) <- var2block$var
  var2block$var <- NULL

  # Compute the aggregated stats for each variable on all bootstrap samples
  res_weights <- res[res$type == "weights", ]
  res_loadings <- res[res$type == "loadings", ]

  stats_weights <- pbapply::pbtapply(
    res_weights$value,
    list(var = res_weights$var, comp = res_weights$comp),
    aggregated_stats, type = "weights"
  )
  stats_loadings <- pbapply::pbtapply(
    res_loadings$value,
    list(var = res_loadings$var, comp = res_loadings$comp),
    aggregated_stats, type = "loadings"
  )
  grid <- expand.grid(
    c(dimnames(stats_weights), list(type = c("weights", "loadings")))
  )

  # Combine the statistics and put them back in the data.frame format
  stats <- rbind(
    matrix(unlist(stats_weights), ncol = 5, byrow = TRUE),
    matrix(unlist(stats_loadings), ncol = 5, byrow = TRUE)
  )
  idx <- which(vapply(stats_loadings, length, 1L) > 0)
  grid <- grid[c(idx, nrow(grid) / 2 + idx), ]

  colnames(stats) <- c("mean", "sd", "pval", "lower_bound", "upper_bound")
  stats <- cbind(stats, grid)
  stats$block <- var2block[stats$var, ]

  ### Compute additional statistics using aggregated statistics
  tail <- qnorm(1 - .05 / 2)
  more_stats <- apply(stats, 1, function(x) {
    if (x["type"] == "weights") {
      ftrans <-
        estimate <- rgcca_res$a[[x["block"]]][x["var"], as.integer(x["comp"])]
    } else {
      estimate <- cor(rgcca_res$blocks[[x["block"]]][, x["var"]],
                      rgcca_res$Y[[x["block"]]][, as.integer(x["comp"])],
                      use = "pairwise.complete.obs"
      )
      ftrans <- 0.5 * log((1 + estimate) / (1 - estimate))
    }
    bootstrap_ratio <- ftrans / as.double(x["sd"])
    return(list(
      estimate = estimate,
      bootstrap_ratio = bootstrap_ratio,
      th_pval = 2 * pnorm(abs(bootstrap_ratio), lower.tail = FALSE),
      th_lower_bound = estimate - as.double(x["sd"]) * tail,
      th_upper_bound = estimate + as.double(x["sd"]) * tail
    ))
  })
  more_stats <- matrix(unlist(more_stats), ncol = 5, byrow = TRUE)
  colnames(more_stats) <- c(
    "estimate", "bootstrap_ratio", "th_pval", "th_lower_bound", "th_upper_bound"
  )
  return(cbind(stats, more_stats))
}
