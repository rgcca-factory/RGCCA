#' Plot fitted rgcca permutation object
#'
#' Plot a fitted rgcca permutation object. The various set of tuning parameters
#' are represented in the x-axis and the the RGCCA objective function - obtained
#' from the both the orginal and permuted blocks - in the y-axis. If type =
#' "zstat" the value of the zstat for the various combinations are reported in
#' the y-axis.
#' @param x A fitted rgcca_permutation object (see
#' \code{\link[RGCCA]{rgcca_permutation}})
#' @param type A character string indicating which criterion to plot
#' (default is 'crit' for the RGCCA criterion or 'zstat' for the pseudo Z-score)
#' @param display_all A boolean indicating is all combinations have to
#' be displayed (default is FALSE).
#' @param show_legend A boolean indicating if the legend is displayed
#' (default is TRUE).
#' @inheritParams plot.rgcca
#' @inheritParams plot_var_2D
#' @inheritParams plot2D
#' @inheritParams get_bootstrap
#' @export
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom stats setNames
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
plot.permutation <- function(x,
                             type = "crit",
                             cex = 1,
                             title = NULL,
                             cex_main = 14 * cex,
                             cex_sub = 12 * cex,
                             cex_point = 3 * cex,
                             cex_lab = 19 * cex,
                             colors = c("red", "grey"),
                             display_all = FALSE,
                             show_legend = TRUE, ...) {
  stopifnot(is(x, "permutation"))
  match.arg(type, c("crit", "zstat"))
  for (i in c("cex", "cex_main", "cex_sub", "cex_point", "cex_lab")) {
    check_integer(i, get(i))
  }
  check_colors(colors)
  if (length(colors) < 2) {
    colors <- rep(colors, 2)
  }
  crit_title <- ifelse(x$call$method %in% c("sgcca", "spls"),
    "SGCCA criterion",
    "RGCCA criterion"
  )

  switch(type,
    "zstat" =  y_title <- "Z-score",
    "crit"  = y_title <- crit_title
  )

  check_ncol(list(x$zstat), 1)

  y <- unlist(x[type])
  N <- nrow(x$penalties)

  df <- setNames(
    data.frame(
      y,
      rep("Other parameter set", N),
      apply(x$penalties, 1, function(x) {
        paste0(round(x, 3), collapse = "/")
      })
    ),
    c(type, "Non_permuted", "label")
  )
  idx_order <- sort(df[[type]], decreasing = FALSE, index.return = TRUE)$ix
  df <- df[idx_order, ]
  best <- which.max(unlist(x["zstat"])[idx_order])

  axis <- function(margin) {
    element_text(
      face = "italic",
      size = cex_lab * 0.75,
      margin = margin
    )
  }

  if (is.null(title)) {
    title <- paste0(
      "Permutation scores (", x$call$n_perms, " runs) \n Best parameters : ",
      df$label[best]
    )
  } else {
    title <- paste0(title, collapse = " ")
  }

  df$label <- factor(df$label, levels = rev(df$label), ordered = TRUE)
  df$Non_permuted[best] <- "Best parameter set"
  breaks <- rev(levels(df$label))
  labels <- as.expression(breaks)
  if (!display_all) {
    jitter <- floor(N / 50)
    breaks[max(best - jitter, 1):min(best + jitter, N)] <- ""
    breaks[best] <- labels[[best]]
  }
  labels[[best]] <- bquote(underline(bold(.(labels[[best]]))))

  p <- ggplot(data = df, mapping = aes_string(x = type, y = "label")) +
    theme_classic() +
    geom_point(aes(shape = Non_permuted, size = Non_permuted)) +
    labs(
      title = title,
      y     = "Combinations",
      x     = y_title
    ) +
    theme_perso(cex, cex_main, cex_sub) +
    theme(
      axis.text = element_text(size = 10, face = "bold"),
      axis.title.x = axis(margin(0, 20, 0, 0)),
      axis.title.y = axis(margin(20, 0, 0, 0)),
      axis.line = element_line(size = 0.5),
      axis.ticks = element_line(size = 0.5),
      axis.ticks.length = unit(2, "mm")
    ) +
    scale_shape_manual("Non Permuted",
      values = c(
        "Best parameter set" = 17,
        "Other parameter set" = 2
      )
    ) +
    scale_size_manual("Non Permuted",
      values = c(
        "Best parameter set" = 4,
        "Other parameter set" = 1
      )
    ) +
    scale_y_discrete(
      labels = labels, breaks = breaks,
      guide = guide_axis(check.overlap = TRUE)
    )
  if (type == "crit") {
    dft <- data.frame(
      combinations = rep(df$label, NCOL(x$permcrit)),
      permcrit = c(x$permcrit[idx_order, ]),
      Permuted = rep(df$Non_permuted, NCOL(x$permcrit))
    )
    p$layers <- c(
      geom_boxplot(
        data = dft,
        aes(
          x = permcrit,
          y = combinations,
          colour = Permuted
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

  p <- p +
    theme(
      plot.margin = margin(5, 0, 0, 0, "mm"),
      legend.position = c(0.7, 0.8)
    )
  if (!show_legend) {
    p <- p + theme(legend.position = "none")
  }
  if (type == "crit") {
    p <- p + expand_limits(x = 0)
  }

  attributes(p)$penalties <- x$penalties[idx_order, ]

  plot(p, ...)
}
