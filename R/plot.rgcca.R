#' Plot for RGCCA
#'
#' Plot different outputs of the results obtained by a rgcca function
#' @param x A fitted RGCCA object (see \code{\link[RGCCA]{rgcca}})
#' @param type A character string: 'sample', 'weight', 'loadings', 'cor_circle',
#' 'both', 'ave' (see details).
#' @param block A numeric vector indicating the blocks to consider.
#' @param comp A numeric vector indicating the components to consider.
#' @param response A vector coloring the points in the "sample" plot.
#' @param title A character string giving the title of the plot.
#' @param cex An integer defining the size of the objects in the plot. Default
#' is one.
#' @param cex_sub An integer defining the font size of the subtitle. Default is
#' 12 * cex.
#' @param cex_main An integer defining the font size of the title. Default is
#' 14 * cex.
#' @param cex_lab An integer defining the font size of the labels. Default is
#' 12 * cex.
#' @param cex_point An integer defining the font size of the points. Default is
#' 3 * cex.
#' @param n_mark An integer defining the maximum number of bars plotted in the
#' "weight" and "loadings" plots.
#' @param colors Colors used in the plots. Default is based on the
#' colorblind-friendly palette recommended in ggplot.
#' @param shapes Shapes used for the points in the "sample" and "cor_circle"
#' plots. Default is the first five shapes used in ggplot.
#' @param ... additional graphical parameters
#' @details
#' \itemize{
#' \item "sample" for sample plot. The blocks (block argument) and components
#' (comp) that will be used on the horizontal and the vertical axes to plot the
#' individuals: (Y[[block[1]]][, comp[1]], Y[[block[2]]][,comp[2]]). Points can
#' be colored according to the response argument. The colors of the points can
#' be modified with the colors argument.
#' \item "weight": barplot of the block weight vector for one
#' specific block/component. The weights are sorted from the highest to
#' the lowest and only the highest are displayed. The number of displayed
#' weights can be set with n_marks.
#' \item "loadings": barplot of the block-loading vector. Variables are sorted
#' in decreasing correlations and only the highest
#' correlations are displayed. The number of displayed correlations can be set
#' with n_marks (defaut value = 30).
#' \item  "cor_circle" for correlation circle.
#' \item "both": displays both sample plot and correlation circle (implemented
#' only for one block and at least when two components are asked (ncomp >= 2).
#' \item "ave": displays the average variance explained for each block.}
#' @examples
#' data(Russett)
#' status <- colnames(Russett)[9:11][apply(Russett[, 9:11], 1, which.max)]
#' X_agric <- as.matrix(Russett[, c("gini", "farm", "rent")])
#' X_ind <- as.matrix(Russett[, c("gnpr", "labo")])
#' X_polit <- as.matrix(Russett[, c("demostab", "dictator")])
#' A <- list(X_agric = X_agric, X_ind = X_ind, X_polit = X_polit)
#' C <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3)
#' fit.rgcca <- rgcca(
#'   blocks = A, connection = C,
#'   tau = rep(1, 3), ncomp = rep(2, 3)
#' )
#'
#' ###############
#' # sample plot #
#' ###############
#' # horizontal axis: First component of the first block
#' # vertical axis: First component of the second block
#' plot(fit.rgcca, type = "sample", block = 1:2, comp = 1, response = status)
#'
#'
#' ######################
#' # all types of plots #
#' ######################
#' # with superblock
#' fit.mcoa <- rgcca(
#'   blocks = A, scheme = "factorial", ncomp = rep(2, 4),
#'   tau = c(1, 1, 1, 0), superblock = TRUE
#' )
#'
#' plot(fit.mcoa, type = "both", response = status)
#' plot(fit.rgcca, type = "loadings")
#' plot(fit.rgcca, type = "weight")
#' plot(fit.rgcca, type = "sample")
#' plot(fit.rgcca, type = "cor_circle")
#' plot(fit.rgcca, type = "ave")
#' @importFrom gridExtra grid.arrange
#' @importFrom ggplot2 ggplot aes
#' @importFrom rlang .data
#' @export
plot.rgcca <- function(x, type = "weight", block = length(x$call$blocks),
                       comp = seq(2),
                       response = as.factor(rep(1, NROW(x$Y[[1]]))),
                       title = NULL, cex = 1, cex_sub = 12 * cex,
                       cex_main = 14 * cex, cex_lab = 12 * cex,
                       cex_point = 3 * cex, n_mark = 30,
                       colors = NULL, shapes = NULL, ...) {
  ### Perform checks and parse arguments
  stopifnot(is(x, "rgcca"))
  type <- tolower(type)
  match.arg(type, c(
    "sample", "cor_circle", "both",
    "ave", "loadings", "weight"
  ))
  lapply(c("cex", "cex_main", "cex_sub", "cex_point", "cex_lab"), function(i) {
    check_integer(i, get(i))
  })
  if (is.null(colors)) {
    colors <- c(
      "#999999", "#E69F00", "#56B4E9", "#009E73",
      "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
    )
  } else {
    check_colors(colors)
  }
  if (is.null(shapes)) {
    shapes <- seq(0, 4)
  } else {
    check_integer("shapes", shapes, min = 0, max = 25, type = "vector")
  }

  # Duplicate comp and block if needed and make sure they take admissible values
  comp <- elongate_arg(comp, seq(2))[seq(2)]
  block <- elongate_arg(block, seq(2))[seq(2)]

  if (type == "ave") comp <- rep(1, 2)
  if (type %in% c("weight", "loadings")) comp <- comp[1]

  lapply(block, function(i) {
    check_blockx("block", i, x$call$blocks)
  })
  Map(
    function(y, z) check_compx(y, y, x$call$ncomp, z), comp, block
  )

  if (missing(response)) {
    if (!is.null(x$call$response)) {
      if (!is.null(x$call$disjunction)) {
        response <- as.factor(x$call$raw[[x$call$response]][, 1])
      } else {
        response <- x$call$blocks[[x$call$response]][, 1]
      }
    }
  } else {
    # Check that response has a good shape
    if ((length(response) != NROW(x$Y[[1]])) || (NCOL(response) > 1)) {
      stop_rgcca(
        "wrong response shape. Response must be a vector of same length as ",
        "the number of subjects seen by RGCCA, i.e., ", NROW(x$Y[[1]])
      )
    }
  }

  # Construct response vector for correlation circle
  if (block[1] > length(x$call$raw)) {
    num_block <- as.factor(unlist(lapply(
      seq_along(x$call$blocks[-block[1]]),
      function(j) rep(names(x$call$blocks)[j], NCOL(x$call$blocks[[j]]))
    )))
  } else {
    num_block <- as.factor(rep(1, NCOL(x$call$blocks[[block[1]]])))
  }

  switch(type,
    # Plot individuals in the projected space
    "sample" = {
      df <- data.frame(
        x$Y[[block[1]]][, comp[1]],
        x$Y[[block[2]]][, comp[2]],
        response = response
      )
      title <- ifelse(missing(title), "Sample space", title)
      plot_function <- plot_sample
    },
    # Plot block columns on a correlation circle
    "cor_circle" = {
      df <- data.frame(cor(
        x$call$blocks[[block[1]]][rownames(x$Y[[block[1]]]), ],
        x$Y[[block[1]]][, comp],
        use = "pairwise.complete.obs"
      ), response = num_block)

      idx <- apply(x$a[[block[1]]][, comp], 1, function(y) {
        any(y != 0)
      })
      df <- df[idx, ]

      title <- ifelse(missing(title), "Correlation circle", title)
      plot_function <- plot_cor_circle
    },
    # Plot side by side individuals in the projected space and block columns on
    # a correlation circle
    "both" = {
      df <- list(
        data.frame(x$Y[[block[1]]][, comp], response = response),
        data.frame(cor(
          x$call$blocks[[block[1]]][rownames(x$Y[[block[1]]]), ],
          x$Y[[block[1]]][, comp],
          use = "pairwise.complete.obs"
        ), response = num_block)
      )

      idx <- apply(x$a[[block[1]]][, comp], 1, function(y) {
        any(y != 0)
      })
      df[[2]] <- df[[2]][idx, ]

      title <- ifelse(
        missing(title), toupper(names(x$call$blocks)[block[1]]), title
      )
      plot_function <- plot_both
    },
    # Plot percentage of AVE per component and per block. If some blocks are not
    # deflated, corrected AVE is reported instead of AVE per component
    "ave" = {
      corrected <- x$call$superblock || !is.null(x$call$response)
      default_title <- ifelse(
        corrected,
        "Corrected Average Variance Explained",
        "Average Variance Explained"
      )
      AVE <- x$AVE$AVE_X_cor
      df <- do.call(rbind, lapply(names(AVE), function(n) {
        data.frame(
          AVE = round(100 * AVE[[n]], digits = 1),
          block = n,
          comp = as.factor(seq(AVE[[n]]))
        )
      }))

      title <- ifelse(missing(title), default_title, title)
      plot_function <- plot_ave
    },
    # Plot the value associated with each individual in the projected space
    "loadings" = {
      df <- data.frame(
        x = cor(
          x$call$blocks[[block[1]]],
          x$Y[[block[1]]][, comp[1]],
          use = "pairwise.complete.obs"
        ),
        y = colnames(x$call$blocks[[block[1]]]),
        response = num_block
      )
      df <- df[order(abs(df$x), decreasing = TRUE), ]
      df <- df[seq(min(n_mark, NROW(df))), ]
      df$y <- factor(df$y, levels = df$y, ordered = TRUE)

      title <- ifelse(missing(title), paste0(
        "Block-loading vector: ",
        names(x$call$blocks)[block[1]], " - comp", comp[1]
      ), title)
      plot_function <- plot_loadings
    },
    # Plot the value associated with each projecting factor
    "weight" = {
      df <- data.frame(
        x = x$a[[block[1]]][, comp[1]],
        y = colnames(x$call$blocks[[block[1]]]),
        response = num_block
      )
      df <- df[order(abs(df$x), decreasing = TRUE), ]
      df <- df[seq(min(n_mark, sum(df$x != 0))), ]
      df$y <- factor(df$y, levels = df$y, ordered = TRUE)

      title <- ifelse(missing(title), paste0(
        "Block-weight vector: ",
        names(x$call$blocks)[block[1]], " - comp", comp[1]
      ), title)
      plot_function <- plot_loadings
    }
  )

  # Construct base theme
  theme_RGCCA <- theme_perso(cex, cex_main, cex_sub, cex_lab)

  # Duplicate colors and shapes to avoid insufficient values in manual scale
  colors <- rep(colors, NROW(df) / length(colors) + 1)
  shapes <- rep(shapes, NROW(df) / length(shapes) + 1)

  # Call plotting function
  p <- plot_function(
    df, title, x, block, comp, theme_RGCCA,
    cex_sub, cex_point, colors, shapes
  )
  if (!is.null(p)) plot(p, ...)
  invisible(p)
}
