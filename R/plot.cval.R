#' Plot cross-validation
#'
#' Plot a cross-validation object (tuning RGCA parameters in 'supervised' mode).
#' The parameters tuned for maximizing RMSE is displayed in the title. In
#' the x-axis, the tuning parameter set. In the y-axis, the average of the
#' cross-validation score. The best parameters are in red by default.
#' @inheritParams plot.rgcca
#' @inheritParams plot.bootstrap
#' @param x A rgcca_cv object (see \link{rgcca_cv})
#' @param type Character string indicating the statistics in the box plots:
#' mean plus or minus standard deviation ("sd", default),
#' \itemize{
#' \item "sd" (default): the middle bar corresponds to the mean and limits of
#' the boxes are given by the mean plus or minus the standard deviation.
#' \item "quantile": the middle bar corresponds to the median and limits of
#' the boxes are given by the 25% and 75% quantiles.
#' }
#' @return A ggplot2 plot object.
#' @examples
#' data("Russett")
#' blocks <- list(
#'   agriculture = Russett[, seq(3)],
#'   industry = Russett[, 4:5],
#'   politic = Russett[, 6:8]
#' )
#' res <- rgcca_cv(blocks,
#'   response = 3, method = "rgcca", par_type = "tau",
#'   par_value = c(0, 0.2, 0.3), n_run = 1, n_cores = 1
#' )
#' plot(res)
#' @export
plot.cval <- function(x, type = "sd",
                      cex = 1,
                      cex_main = 14 * cex,
                      cex_sub = 12 * cex,
                      cex_point = 3 * cex,
                      cex_lab = 12 * cex,
                      display_order = TRUE, ...) {
  ### Perform checks and parse params
  stopifnot(is(x, "cval"))
  type <- match.arg(type, c("quantile", "sd"))

  ### Build data frame
  df <- data.frame(
    combinations = x$stats$combinations,
    ymin = apply(x$cv, 1, min),
    ymax = apply(x$cv, 1, max)
  )
  if (type == "quantile") {
    df$middle <- x$stats$median
    df$lower <- x$stats$Q1
    df$upper <- x$stats$Q3
  } else {
    df$middle <- x$stats$mean
    df$lower <- x$stats$mean - x$stats$sd
    df$upper <- x$stats$mean + x$stats$sd
  }

  best <- which(apply(
    x$penalties, 1, function(z) identical(z, x$bestpenalties)
  ))
  df$category <- rep("param", NROW(x$cv))
  df$category[best] <- "best_param"

  labels <- as.expression(df$combinations)
  labels[[best]] <- bquote(underline(bold(.(labels[[best]]))))

  if (display_order) {
    idx_order <- sort(df$middle, decreasing = FALSE, index.return = TRUE)$ix
    df <- df[idx_order, ]
  }
  df$combinations <- factor(
    df$combinations, levels = df$combinations, ordered = TRUE
  )

  df_points <- data.frame(
    combinations = rep(df$combinations, NCOL(x$cv)),
    y = c(x$cv),
    category = rep(df$category, NCOL(x$cv))
  )

  ### Prepare plot
  n_run_str <- ifelse(
    x$n_run > 1, paste0(" and ", x$n_run, " run(s)"), ""
  )
  validation_str <- ifelse(
    x$validation == "kfold",
    paste0("kfold: with ", x$k, " folds", n_run_str, ""),
    "leave-one-out"
  )

  title <- paste0(
    "Cross-validated ", x$metric, " (", validation_str, ")\nBest parameters: ",
    df$combinations[best]
  )
  xlab <- paste0("Tuning parameter sets (", x$par_type, ")")
  ylab <- paste("Mean", x$metric)

  ### Construct plot
  p <- ggplot(data = df, aes(
    x = .data$combinations,
    group = .data$combinations,
    color = .data$category
  )) +
    ggplot2::coord_flip() +
    ggplot2::scale_color_manual(
      values = c("param" = "grey", "best_param" = "black")
    ) +
    ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab) +
    ggplot2::scale_x_discrete(
      labels = labels, breaks = df$combinations,
      guide = ggplot2::guide_axis(check.overlap = TRUE)
    ) +
    ggplot2::geom_boxplot(
      aes(
        ymin = .data$ymin, lower = .data$lower, middle = .data$middle,
        upper = .data$upper, ymax = .data$ymax
      ),
      stat = "identity"
    ) +
    ggplot2::geom_jitter(
      data = df_points, aes(
        x = .data$combinations, y = .data$y, color = .data$category
      ), size = .5 * cex_point,
      position = ggplot2::position_jitter(height = 0.01, width = 0.2)
    ) +
    ggplot2::ggtitle(title) +
    theme_perso(cex, cex_main, cex_sub, cex_lab) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = .8 * cex_sub),
      legend.position = "none"
    )
  plot(p, ...)
  invisible(p)
}
