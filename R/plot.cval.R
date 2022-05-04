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
#' \item "stderr": the middle bar corresponds to the mean and limits of
#' the boxes are given by the mean plus or minus the standard deviation divided
#' by the square roots of the number of folds.
#' \item "quantile": the middle bar corresponds to the median and limits of
#' the boxes are given by the 5% and 95% quantiles.
#' \item "points": box plots are removed and only the points are kept.}
#' @param colors Colors used in the plots. Default is black and red.
#' @examples
#' data("Russett")
#' blocks <- list(
#'   agriculture = Russett[, seq(3)],
#'   industry = Russett[, 4:5],
#'   politic = Russett[, 6:11]
#' )
#' res <- rgcca_cv(blocks,
#'   response = 3, method = "rgcca", par_type = "tau",
#'   par_value = c(0, 0.2, 0.3), n_run = 1, n_cores = 1
#' )
#' plot(res)
#' @export
plot.cval <- function(x, type = "sd",
                      cex = 1, cex_main = 14 * cex,
                      cex_sub = 10 * cex, cex_lab = 10 * cex,
                      colors = c("grey", "red"),
                      display_order = TRUE, ...) {
  ### Perform checks and parse params
  stopifnot(is(x, "cval"))
  match.arg(type, c("quantile", "sd", "stderr", "points"))
  for (i in c("cex", "cex_main", "cex_sub", "cex_lab")) {
    check_integer(i, get(i))
  }
  check_colors(colors)
  colors <- elongate_arg(colors, seq(2))

  ### Build data frame
  ymin <- apply(x$cv, 1, min)
  ymax <- apply(x$cv, 1, max)
  ymean <- apply(x$cv, 1, mean)
  ymed <- apply(x$cv, 1, median)

  switch(type,
    "quantile" = {
      middle <- ymed
      lower <- apply(x$cv, 1, quantile, 0.05)
      upper <- apply(x$cv, 1, quantile, 0.95)
    },
    "sd" = {
      middle <- ymean
      lower <- middle - apply(x$cv, 1, sd)
      upper <- middle + apply(x$cv, 1, sd)
    },
    "stderr" = {
      middle <- ymean
      lower <- middle - apply(x$cv, 1, function(y) {
        sd(y) / sqrt(length(y))
      })
      upper <- middle + apply(x$cv, 1, function(y) {
        sd(y) / sqrt(length(y))
      })
    },
    "points" = {
      middle <- ymean
      lower <- apply(x$cv, 1, min)
      upper <- apply(x$cv, 1, max)
    }
  )

  if (length(x$call$blocks) > 3) {
    combinations <- paste("Set ", seq(NROW(x$penalties)))
  } else {
    combinations <- apply(
      format(x$penalties, digits = 2), 1, paste0,
      collapse = "/"
    )
  }

  if (display_order) {
    idx_order <- sort(middle, decreasing = FALSE, index.return = TRUE)$ix
    combinations <- factor(
      combinations,
      levels = combinations[idx_order], ordered = TRUE
    )
  } else {
    combinations <- factor(combinations, levels = combinations, ordered = TRUE)
  }

  category <- rep("param", NROW(x$cv))
  category[which(apply(
    x$penalties, 1, function(z) identical(z, x$bestpenalties)
  ))] <- "best_param"

  df <- data.frame(
    combinations = combinations,
    ymin = ymin,
    lower = lower,
    middle = middle,
    upper = upper,
    ymax = ymax,
    category = category
  )

  df_points <- data.frame(
    combinations = rep(combinations, NCOL(x$cv)),
    y = c(x$cv),
    category = rep(category, NCOL(x$cv))
  )

  ### Prepare plot
  n_run_str <- ifelse(
    x$call$n_run > 1, paste0(" and ", x$call$n_run, " run(s)"), ""
  )
  validation_str <- ifelse(
    x$call$validation == "kfold",
    paste0("kfold: with ", x$call$k, " folds", n_run_str, ""),
    "leave-one-out"
  )

  title <- paste0(
    "Cross-validated error (", validation_str, ")\nBest parameters: ",
    paste(round(x$bestpenalties, digits = 2), collapse = "/")
  )
  xlab <- paste0("Tuning parameter sets (", x$call$par_type, ")")
  ylab <- "Mean error"

  ### Construct plot
  p <- ggplot(data = df, aes(
    x = .data$combinations,
    group = .data$combinations,
    fill = .data$category
  )) +
    ggplot2::geom_point(data = df_points, aes(
      x = .data$combinations, y = .data$y, color = .data$category
    )) +
    ggplot2::scale_fill_manual(
      values = c("param" = colors[1], "best_param" = colors[2])
    ) +
    ggplot2::scale_color_manual(
      values = c("param" = colors[1], "best_param" = colors[2])
    ) +
    ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab)
  if (type != "points") {
    p <- p +
      ggplot2::geom_boxplot(
        aes(
          ymin = .data$ymin, lower = .data$lower, middle = .data$middle,
          upper = .data$upper, ymax = .data$ymax
        ),
        stat = "identity"
      )
  }

  # Set theme
  p <- p + ggplot2::ggtitle(title) + theme_perso(cex, cex_main, cex_sub) +
    ggplot2::theme(
      axis.text = ggplot2::element_text(size = 10, face = "bold"),
      axis.line = ggplot2::element_line(size = 0.5),
      axis.ticks = ggplot2::element_line(size = 0.5),
      axis.ticks.length = ggplot2::unit(2, "mm"),
      legend.position = "none",
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
    )
  plot(p, ...)
}
