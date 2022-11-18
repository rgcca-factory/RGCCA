#' Plot a fitted bootstrap object
#'
#' Plot the results of a fitted bootstrap object. Each block variable is shown
#' along with its associated bootstrap confidence interval and stars reflecting
#' the p-value of assigning a strictly positive or negative weight to this
#' block variable.
#' @inheritParams plot.rgcca
#' @param x A fitted bootstrap object (see \code{\link[RGCCA]{bootstrap}})
#' @param type Character string indicating the bootstrapped object to plot:
#' block-weight vectors ("weights", default) or block-loading vectors
#' ("loadings").
#' @param empirical A logical value indicating if the bootstrap confidence
#' intervals and p-values are derived from the empirical distribution.
#' (defaut: TRUE)
#' @param n_mark An integer defining the maximum number of variables plotted.
#' @param show_sign A logical for showing significance levels.
#'
#' @examples
#' data("Russett")
#' blocks <- list(
#'   agriculture = Russett[, seq(3)],
#'   industry = Russett[, 4:5],
#'   politic = Russett[, 6:11]
#' )
#' fit.rgcca <- rgcca(blocks, ncomp = 2, method = "rgcca", tau = 1)
#' fit.boot <- bootstrap(fit.rgcca, n_boot = 20, n_cores = 1)
#' plot(fit.boot, type = "weight", block = 1, comp = 1)
#' @export
plot.bootstrap <- function(x, block = seq_along(x$rgcca$call$raw),
                           comp = 1, type = "weights",
                           empirical = TRUE, n_mark = 30,
                           display_order = TRUE,
                           show_sign = TRUE, title = NULL,
                           cex = 1, cex_sub = 12 * cex,
                           cex_main = 14 * cex, cex_lab = 12 * cex,
                           cex_point = 3 * cex, colors = NULL, ...) {
  ### Perform checks and parse arguments
  stopifnot(is(x, "bootstrap"))
  type <- match.arg(type, c("weights", "loadings"))
  lapply(block, function(i) check_blockx("block", i, x$rgcca$call$raw))
  check_integer("n_mark", n_mark)

  if (is.null(colors)) {
    colors <- c(
      "#999999", "#E69F00", "#56B4E9", "#009E73",
      "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
    )
  } else {
    check_colors(colors)
  }

  ### Build data frame
  df <- Reduce(rbind, lapply(
    block,
    function(b) {
      get_bootstrap(
        b = x, type = type,
        block = b,
        comp = comp,
        empirical = empirical
      )
    }
  ))
  df$response <- as.factor(unlist(lapply(block, function(j) {
    rep(names(x$rgcca$call$blocks)[j], NCOL(x$rgcca$call$blocks[[j]]))
  })))

  if (display_order) {
    df <- df[order(abs(df$estimate), decreasing = TRUE), ]
  }
  df <- df[df[, "sd"] != 0, ]

  n_mark <- min(n_mark, NROW(df))
  df <- data.frame(df, order = seq(NROW(df), 1))[seq(n_mark), ]

  significance <- rep("", n_mark)
  significance[df$pval < 1e-3] <- "***"
  significance[df$pval >= 1e-3 & df$pval < 1e-2] <- "**"
  significance[df$pval >= 1e-2 & df$pval < 5e-2] <- "*"

  df$sign <- factor(significance,
    levels = c(labels = c("***", "**", "*", ""))
  )

  ### Prepare plot
  block_name <- ifelse(
    length(block) > 1,
    "",
    paste0("(", names(x$rgcca$call$blocks)[block], ")")
  )

  title <- ifelse(is.null(title),
    paste0(
      "Bootstrap confidence interval ",
      block_name, "\n (",
      type, " - ",
      ncol(x$bootstrap[[1]][[1]][[1]]),
      " bootstrap samples - comp ", comp, ")"
    ),
    title
  )

  ### Construct plot
  if (length(block) > 1) {
    p <- ggplot(
      df,
      aes(x = .data$order, y = .data$estimate, color = .data$response)
    ) +
      ggplot2::scale_color_manual(values = colors) +
      ggplot2::labs(color = "Block")
  } else {
    p <- ggplot(df, aes(x = .data$order, y = .data$estimate))
  }

  p <- p +
    ggplot2::geom_point(size = .6 * cex_point) +
    ggplot2::geom_errorbar(aes(
      ymin = .data$lower_bound,
      ymax = .data$upper_bound,
      width = .2 * cex_point
    ), size = .2 * cex_point) +
    ggplot2::coord_flip() +
    theme_perso(cex, cex_main, cex_sub, cex_lab) +
    ggplot2::labs(title = title, x = "", y = "") +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(
        size = cex_sub,
        face = "italic",
        color = "gray40"
      ),
      axis.text.x = ggplot2::element_text(
        size = cex_sub,
        face = "italic",
        color = "gray40"
      ),
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(5, 0, 0, 0, "mm")
    ) +
    ggplot2::scale_x_continuous(
      breaks = df$order,
      labels = rownames(df)
    ) +
    ggplot2::geom_hline(yintercept = 0, lty = "longdash")
  if (show_sign) {
    p <- p + ggplot2::geom_text(
      aes(label = .data$sign),
      nudge_x = 0.1, size = 2 * cex_point, show.legend = FALSE
    )
  }

  plot(p, ...)
}
