#' Plot a fitted bootstrap object
#'
#' Plot the results of a fitted bootstrap object. A bar plot is constructed
#' where each block variable is shown along with its associated bootstrap
#' confidence interval and a color reflecting the p-value of assigning a
#' strictly positive or negative weight to this block variable.
#' @inheritParams plot.rgcca
#' @param x A fitted bootstrap object (see \code{\link[RGCCA]{bootstrap}})
#' @param type Character string indicating the bootstrapped object to plot:
#' block-weight vectors ("weight", default) or block-loading vectors
#' ("loadings").
#' @param empirical A logical value indicating if the bootstrap confidence
#' intervals and p-values are derived from the empirical distribution.
#' (defaut: TRUE)
#' @param display_order A logical value for ordering the variables.
#' @param n_mark An integer defining the maximum number of bars plotted.
#' @param colors Colors used in the plots. Default is a grey scaled palette.
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
#' @importFrom grDevices grey.colors
#' @export
plot.bootstrap <- function(x, block = length(x$rgcca$call$blocks),
                           comp = 1, type = "weight",
                           empirical = TRUE, n_mark = 30,
                           display_order = TRUE,
                           colors = grey.colors(6)[2:6], title = NULL,
                           cex = 1, cex_sub = 12 * cex,
                           cex_main = 14 * cex, cex_lab = 12 * cex,
                           cex_point = 3 * cex, ...) {
  ### Perform checks and parse arguments
  stopifnot(is(x, "bootstrap"))
  check_blockx("block", block, x$rgcca$call$blocks)
  check_integer("n_mark", n_mark)
  check_colors(colors)
  colors <- elongate_arg(colors, seq(5))[seq(5)]

  ### Build data frame
  df <- get_bootstrap(
    x,
    type = type,
    comp,
    block = block,
    empirical = empirical,
    display_order = display_order
  )
  df <- df[df[, "sd"] != 0, ]

  n_mark <- min(n_mark, NROW(df))
  df <- data.frame(df, order = seq(NROW(df), 1))[seq(n_mark), ]

  significance <- rep("", n_mark)
  significance[df$pval < 1e-3] <- "< 0.001"
  significance[df$pval >= 1e-3 & df$pval < 1e-2] <- "< 0.01"
  significance[df$pval >= 1e-2 & df$pval < 5e-2] <- "< 0.05"
  significance[df$pval >= 5e-2 & df$pval < 1e-1] <- "< 0.1"
  significance[df$pval >= 1e-1] <- "> 0.1"

  df$sign <- factor(significance,
    levels = c(labels = c(
      "< 0.001", "< 0.01",
      "< 0.05", "< 0.1", "> 0.1"
    ))
  )

  ### Prepare plot
  title <- ifelse(is.null(title),
    paste0(
      "Bootstrap confidence interval (",
      names(x$rgcca$call$blocks)[block], ")\n (",
      type, " - ",
      ncol(x$bootstrap[[1]][[1]][[1]]),
      " bootstrap samples - comp ", comp, ")"
    ),
    title
  )

  ### Construct plot
  p <- ggplot(
    df,
    aes(x = .data$order, y = .data$estimate, fill = .data$sign)
  ) +
    ggplot2::geom_bar(stat = "identity") +
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
    ggplot2::scale_fill_manual(
      values = colors,
      labels = c(
        "< 0.001", "< 0.01", "< 0.05",
        "< 0.1", "> 0.1"
      ),
      drop = FALSE,
      name = "Signif."
    )

  if (n_mark <= 50) {
    p <- p + ggplot2::geom_errorbar(aes(
      ymin = .data$lower_bound,
      ymax = .data$upper_bound,
      width = 0.5
    ))
  }

  plot(p, ...)
}
