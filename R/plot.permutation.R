#' Plot a fitted rgcca permutation object
#'
#' Plot a fitted rgcca permutation object. The various set of tuning parameters
#' are represented in the x-axis and the RGCCA objective function - obtained
#' from the both the orginal and permuted blocks - in the y-axis. If type =
#' "zstat" the value of the zstat for the various parameter sets are reported in
#' the y-axis.
#' @inheritParams plot.rgcca
#' @param x A fitted rgcca_permutation object (see
#' \code{\link[RGCCA]{rgcca_permutation}}).
#' @param type A character string indicating which criterion to plot (default
#' is 'crit' for the RGCCA criterion or 'zstat' for the pseudo Z-score).
#' @param display_all A boolean indicating if all parameter sets have to
#' be displayed (default is FALSE).
#' @param show_legend A boolean indicating if legend should
#' be shown (default is TRUE).
#' @examples
#' data(Russett)
#' A <- list(
#'   agriculture = Russett[, seq(3)],
#'   industry = Russett[, 4:5],
#'   politic = Russett[, 6:11]
#' )
#'
#' perm.out <- rgcca_permutation(A, par_type = "tau", n_perms = 2, n_cores = 1)
#' print(perm.out)
#' plot(perm.out)
#'
#' perm.out <- rgcca_permutation(A,
#'   par_type = "sparsity",
#'   n_perms = 5, n_cores = 1
#' )
#' print(perm.out)
#' plot(perm.out, type = "zstat")
#' @export
plot.permutation <- function(x,
                             type = "crit",
                             cex = 1,
                             title = NULL,
                             cex_main = 14 * cex,
                             cex_sub = 12 * cex,
                             cex_point = 3 * cex,
                             cex_lab = 19 * cex,
                             display_all = FALSE,
                             show_legend = FALSE, ...) {
  ### Perform checks and parse params
  stopifnot(is(x, "permutation"))
  match.arg(type, c("crit", "zstat"))
  for (i in c("cex", "cex_main", "cex_sub", "cex_point", "cex_lab")) {
    check_integer(i, get(i))
  }

  ### Build data frame
  if (length(x$call$blocks) > 5) {
    combinations <- paste("Set ", seq_along(x$pvals))
  } else {
    combinations <- apply(
      format(x$penalties, digits = 2), 1, paste0,
      collapse = "/"
    )
  }

  df <- data.frame(
    x = x[[type]],
    label = "Other parameter set",
    combinations = combinations
  )

  # Reorder dataframe according to the quantity of interest
  idx_order <- sort(df$x, decreasing = FALSE, index.return = TRUE)$ix
  df <- df[idx_order, ]

  # Mark the best parameter set
  best <- which(apply(
    x$penalties[idx_order, ], 1, function(z) identical(z, x$bestpenalties)
  ))
  df$label[best] <- "Best parameter set"

  df$combinations <- factor(
    df$combinations,
    levels = rev(df$combinations), ordered = TRUE
  )

  ### Prepare plot
  crit_title <- ifelse(x$call$method %in% c("sgcca", "spls"),
    "SGCCA criterion",
    "RGCCA criterion"
  )
  xlab <- ifelse(type == "zstat", "Z-score", crit_title)
  ylab <- paste0("Tuning parameter sets (", x$call$par_type, ")")

  breaks <- rev(levels(df$combinations))
  labels <- as.expression(breaks)
  if (!display_all) {
    N <- nrow(x$penalties)
    jitter <- floor(N / 50)
    breaks[max(best - jitter, 1):min(best + jitter, N)] <- ""
    breaks[best] <- labels[[best]]
  }
  labels[[best]] <- bquote(underline(bold(.(labels[[best]]))))

  title <- ifelse(
    missing(title),
    paste0(
      "Permutation scores (", x$call$n_perms, " runs) \n Best parameters: ",
      df$combinations[best]
    ),
    title
  )

  ### Construct plot
  # Main plot (values of the quantity of interest per set of parameters)
  p <- ggplot(
    data = df, mapping = aes_(x = quote(x), y = quote(combinations))
  ) +
    theme_classic() +
    geom_point(aes_(shape = quote(label), size = quote(label))) +
    labs(
      title = title,
      y     = ylab,
      x     = xlab
    ) +
    theme_perso(cex, cex_main, cex_sub) +
    theme(
      axis.text = element_text(size = 10, face = "bold"),
      axis.line = element_line(size = 0.5),
      axis.ticks = element_line(size = 0.5),
      axis.ticks.length = unit(2, "mm")
    ) +
    scale_shape_manual("Non permuted",
      values = c(
        "Best parameter set" = 17,
        "Other parameter set" = 2
      )
    ) +
    scale_size_manual("Non permuted",
      values = c(
        "Best parameter set" = 4,
        "Other parameter set" = 1
      )
    ) +
    scale_y_discrete(
      labels = labels, breaks = breaks,
      guide = guide_axis(check.overlap = TRUE)
    )

  # Second part of the plot (boxplots of the permuted criteria)
  if (type == "crit") {
    df2 <- data.frame(
      combinations = rep(df$combinations, NCOL(x$permcrit)),
      x = c(x$permcrit[idx_order, ]),
      label = rep(df$label, NCOL(x$permcrit))
    )
    p$layers <- c(
      geom_boxplot(
        data = df2,
        aes_(
          x = quote(x),
          y = quote(combinations),
          colour = quote(label)
        ),
        size = 0.8
      ),
      p$layers
    )
    p <- p + scale_colour_manual("Permuted", values = c(
      "Best parameter set" = "black",
      "Other parameter set" = "grey"
    ))
  }

  # Move legend
  if (!show_legend) {
    p <- p + theme(legend.position = "none")
  } else {
    p <- p +
      theme(
        legend.position = c(.9, 1.),
        legend.justification = c("right", "top")
      )
  }

  plot(p, ...)
}
