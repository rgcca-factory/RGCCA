#' Plot for RGCCA
#'
#' Plot different outputs of the results obtained by a rgcca function
#' @param x A fitted RGCCA object (see \code{\link[RGCCA]{rgcca}})
#' @param type A character string: 'samples', 'weights', 'loadings',
#' 'cor_circle', both', 'ave' (see details).
#' @param block A numeric corresponding to the block(s) to plot.
#' @param comp A numeric vector indicating the components to consider.
#' @param response A vector coloring the points in the "samples" plot.
#' @param display_order A logical value for ordering the variables.
#' @param title A character string giving the title of the plot.
#' @param cex A numeric defining the size of the objects in the plot. Default
#' is one.
#' @param cex_sub A numeric defining the font size of the subtitle. Default is
#' 12 * cex.
#' @param cex_main A numeric defining the font size of the title. Default is
#' 14 * cex.
#' @param cex_lab A numeric defining the font size of the labels. Default is
#' 12 * cex.
#' @param cex_point A numeric defining the font size of the points. Default is
#' 3 * cex.
#' @param n_mark An integer defining the maximum number of bars plotted in the
#' "weights" and "loadings" plots.
#' @param sample_colors Colors used to color samples.
#' @param sample_shapes Shapes used for the sample points.
#' @param var_colors Colors used to color variable weights or correlations
#' with canonical components.
#' @param var_shapes Shapes used for the points associated to variable weights
#' or correlations with canonical components.
#' @param AVE_colors Colors used in the AVE plot.
#' @param show_labels A logical value for showing the labels in plots "samples"
#' and "cor_circle".
#' @param repel A logical value for repelling text labels from each other.
#' @param display_blocks A numeric corresponding to the block(s) to display in
#' the correlation_circle.
#' @param expand A numeric that scales the weights associated to the block
#' variables in the biplot.
#' @param show_arrows A logical, if TRUE, arrows are shown in the biplot.
#' Default is FALSE.
#' @param ... additional graphical parameters
#' @details
#' \itemize{
#' \item "samples" for sample plot. The blocks (block argument) and components
#' (comp) that will be used on the horizontal and the vertical axes to plot the
#' individuals: (Y[[block[1]]][, comp[1]], Y[[block[2]]][,comp[2]]). Points can
#' be colored according to the response argument. The colors of the points can
#' be modified with the colors argument.
#' \item "weights": barplot of the block weight vector for one
#' specific block/component. The weights are sorted from the highest to
#' the lowest and only the highest are displayed. The number of displayed
#' weights can be set with n_marks.
#' \item "loadings": barplot of the block-loading vector. Variables are sorted
#' in decreasing correlations and only the highest
#' correlations are displayed. The number of displayed correlations can be set
#' with n_marks (defaut value = 30).
#' \item  "cor_circle" for correlation circle. It represents the correlation
#' between the component corresponding to the first element of the block
#' argument, and the variables of the block corresponding to the blocks
#' specified by the argument display_blocks.
#' \item "both": displays both sample plot and correlation circle (implemented
#' only for one block and at least when two components are asked (ncomp >= 2).
#' \item "biplot": displays on the same plot the samples and the variables
#' used to compute the components.
#' \item "ave": displays the average variance explained for each block.}
#' @return A ggplot2 plot object.
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
#' plot(fit.rgcca, type = "samples", block = 1:2, comp = 1, response = status)
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
#' plot(fit.rgcca, type = "biplot")
#' plot(fit.rgcca, type = "ave")
#' @importFrom gridExtra grid.arrange
#' @importFrom ggplot2 ggplot aes
#' @importFrom ggrepel geom_text_repel
#' @importFrom rlang .data
#' @export
plot.rgcca <- function(x, type = "weights", block = seq_along(x$call$blocks),
                       comp = c(1, 2),
                       response = as.factor(rep(1, NROW(x$Y[[1]]))),
                       display_order = TRUE,
                       title = NULL, cex = 1, cex_sub = 12 * cex,
                       cex_main = 14 * cex, cex_lab = 12 * cex,
                       cex_point = 3 * cex, n_mark = 30,
                       sample_colors = NULL, sample_shapes = NULL,
                       var_colors = NULL, var_shapes = NULL,
                       AVE_colors = NULL,
                       show_labels = TRUE, repel = FALSE,
                       display_blocks = seq_along(x$call$blocks),
                       expand = 1, show_arrows = TRUE, ...) {
  ### Define data.frame generating functions
  df_sample <- function(x, block, comp, response, obj = "Y") {
    data.frame(
      x[[obj]][[block[1]]][, comp[1]],
      x[[obj]][[block[2]]][, comp[2]],
      response = response
    )
  }

  df_cor <- function(x, block, comp, num_block, display_blocks) {
    df <- data.frame(
      x = do.call(rbind, Map(function(i, j) {
        cor(
          x$blocks[[i]][rownames(x$Y[[j]]), ],
          x$Y[[j]][, comp],
          use = "pairwise.complete.obs"
        )
      }, display_blocks, block)),
      response = num_block,
      y = do.call(c, lapply(x$blocks[display_blocks], colnames))
    )

    idx <- apply(do.call(
      rbind,
      lapply(x$a[display_blocks], function(z) z[, comp, drop = FALSE])
    ), 1, function(y) any(y != 0))
    return(df[idx, ])
  }

  df_weight <- function(x, block, comp, num_block, display_order) {
    df <- data.frame(
      x = unlist(lapply(x$a[block], function(z) z[, comp[1]])),
      y = do.call(c, lapply(x$blocks[block], colnames)),
      response = num_block
    )
    df <- df[df$x != 0, ]
    if (display_order) {
      df <- df[order(abs(df$x), decreasing = TRUE), ]
    }
    df <- df[seq(min(n_mark, nrow(df))), ]
    df$y <- factor(df$y, levels = df$y, ordered = TRUE)
    return(df)
  }

  df_AVE <- function(x) {
    AVE <- x$AVE$AVE_X_cor
    df <- do.call(rbind, lapply(names(AVE), function(n) {
      data.frame(
        AVE = round(100 * AVE[[n]], digits = 1),
        block = n,
        comp = as.factor(seq(AVE[[n]]))
      )
    }))
    df$label <- df$AVE
    df$label[df$AVE == 0] <- ""
    return(df)
  }

  ### Perform checks and parse arguments
  stopifnot(is(x, "rgcca"))
  check_boolean(display_order)
  check_boolean(show_labels)
  check_boolean(show_arrows)
  check_integer("expand", expand, float = TRUE, min = -Inf)
  check_integer("comp", comp, type = "vector")
  type <- tolower(type)
  type <- match.arg(type, c(
    "samples", "cor_circle", "both",
    "ave", "loadings", "weights", "biplot"
  ))
  sample_colors <- check_colors(sample_colors, type = "samples")
  var_colors <- check_colors(var_colors, type = "variables")
  AVE_colors <- check_colors(AVE_colors, type = "AVE")
  if (is.null(sample_shapes)) {
    sample_shapes <- seq(0, 4)
  } else {
    check_integer(
      "sample_shapes", sample_shapes, min = 0, max = 25, type = "vector"
    )
  }
  if (is.null(var_shapes)) {
    var_shapes <- seq(15, 19)
  } else {
    check_integer("var_shapes", var_shapes, min = 0, max = 25, type = "vector")
  }

  # Duplicate comp and block if needed and make sure they take admissible values
  comp <- elongate_arg(comp, seq(2))[seq(2)]
  block <- elongate_arg(block, seq(2))

  lapply(display_blocks, function(i) {
    check_blockx("display_blocks", i, x$call$blocks)
  })

  response_block <- FALSE

  switch(type,
    "samples" = {
      block <- block[seq(2)]
      display_blocks <- block
    },
    "cor_circle" = {
      block <- ifelse(x$call$superblock, length(x$call$blocks) + 1, block[1])
      response_block <- !is.null(x$call$response) && (block == x$call$response)
    },
    "both" = {
      block <- rep(ifelse(
        x$call$superblock, length(x$call$blocks) + 1, block[1]
      ), 2)
      response_block <- !is.null(x$call$response) &&
        (block[1] == x$call$response)
    },
    "ave" = {
      comp <- 1
    },
    "loadings" = {
      block <- unique(block)
      comp <- comp[1]
      display_blocks <- block
    },
    "weights" = {
      block <- unique(block)
      comp <- comp[1]
      display_blocks <- block
    },
    "biplot" = {
      block <- rep(ifelse(
        x$call$superblock, length(x$call$blocks) + 1, block[1]
      ), 2)
      display_blocks <- block[1]
      response_block <- !is.null(x$call$response) &&
        (block[1] == x$call$response)
    }
  )

  if (response_block) {
    stop_rgcca(
      "Drawing a correlation circle or a biplot with block = ", block ,
      " is not allowed, because the response components are not orthogonal."
    )
  }

  lapply(block, function(i) {
    check_blockx("block", i, x$blocks)
  })
  comp <- unlist(Map(
    function(y, z) min(y, x$call$ncomp[z]), comp, block
  ))[seq_along(comp)]

  if (length(unique(comp)) == 1 & type %in% c("both", "cor_circle", "biplot")) {
    stop_rgcca(
      "Drawing a correlation circle or a biplot with a single component ",
      "is not allowed, please provide different values for comp."
    )
  }

  if (missing(response)) {
    if (!is.null(x$call$response)) {
      if (x$opt$disjunction) {
        response <- as.factor(x$call$blocks[[x$call$response]][, 1])
      } else {
        response <- x$blocks[[x$call$response]][, 1]
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

  # Construct response vector for correlation circle, weights and loadings
  if (block[1] == length(x$call$blocks) + 1) {
    num_block <- as.factor(unlist(lapply(
      seq_along(x$blocks[-block[1]]),
      function(j) rep(names(x$blocks)[j], NCOL(x$blocks[[j]]))
    )))
  } else {
    num_block <- as.factor(unlist(lapply(
      display_blocks,
      function(j) rep(names(x$blocks)[j], NCOL(x$blocks[[j]]))
    )))
  }

  switch(type,
    # Plot individuals in the projected space
    "samples" = {
      df <- df_sample(x, block, comp, response)
      title <- ifelse(missing(title), "Sample space", title)
      plot_function <- plot_sample
    },
    # Plot block columns on a correlation circle
    "cor_circle" = {
      df <- df_cor(x, block[1], comp, num_block, display_blocks)

      title <- ifelse(missing(title), "Correlation circle", title)
      plot_function <- plot_cor_circle
    },
    # Plot side by side individuals in the projected space and block columns on
    # a correlation circle
    "both" = {
      df <- list(
        df_sample(x, block, comp, response),
        df_cor(x, block[1], comp, num_block, display_blocks)
      )

      title <- ifelse(missing(title), "", title)
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
      df <- df_AVE(x)

      title <- ifelse(missing(title), default_title, title)
      plot_function <- plot_ave
    },
    # Plot the value associated with each individual in the projected space
    "loadings" = {
      df <- df_cor(x, block, comp, num_block, display_blocks)
      if (display_order) {
        df <- df[order(abs(df[, 1]), decreasing = TRUE), ]
      }
      df <- df[seq(min(n_mark, NROW(df))), ]
      df$y <- factor(df$y, levels = df$y, ordered = TRUE)

      block_name <- ifelse(
        length(unique(block)) == 1,
        paste(":", names(x$blocks)[block[1]]),
        ""
      )
      title <- ifelse(missing(title), paste0(
        "Block-loading vector", block_name, " - comp", comp[1]
      ), title)
      plot_function <- plot_loadings
    },
    # Plot the value associated with each projecting factor
    "weights" = {
      df <- df_weight(x, block, comp, num_block, display_order)

      block_name <- ifelse(
        length(unique(block)) == 1,
        paste(":", names(x$blocks)[block[1]]),
        ""
      )
      title <- ifelse(missing(title), paste0(
        "Block-weight vector", block_name, " - comp", comp[1]
      ), title)
      plot_function <- plot_loadings
    },
    # Plot individuals in the projected space and variables
    # used for the projection
    "biplot" = {
      df <- list(
        Y = df_sample(x, block, comp, response, "Y"),
        a = df_sample(x, block, comp, num_block, "a")
      )

      # Rescale weigths
      var_tot <- sum(diag(var(x$blocks[[block[1]]], na.rm = TRUE)))
      a <- data.matrix(df$a[, c(1, 2)]) %*% diag(sqrt(
        var_tot * x$AVE$AVE_X[[block[1]]][comp]
      ))

      ratio <- min(
        (max(df$Y[, 1]) - min(df$Y[, 1])) / (max(a[, 1]) - min(a[, 1])),
        (max(df$Y[, 2]) - min(df$Y[, 2])) / (max(a[, 2]) - min(a[, 2]))
      ) * expand

      df$a[, c(1, 2)] <- a * ratio
      df$a <- df$a[df$a[, 1] != 0 | df$a[, 2] != 0, ]

      title <- ifelse(missing(title), "Biplot", title)
      plot_function <- plot_biplot
    }
  )

  # Construct base theme
  theme_RGCCA <- theme_perso(cex, cex_main, cex_sub, cex_lab)

  # Duplicate colors and shapes to avoid insufficient values in manual scale
  sample_colors <- rep(
    sample_colors,
    ifelse(is.list(df), NROW(df[[1]]), NROW(df)) / length(sample_colors) + 1
  )
  sample_shapes <- rep(
    sample_shapes,
    ifelse(is.list(df), NROW(df[[1]]), NROW(df)) / length(sample_shapes) + 1
  )
  var_colors <- rep(
    var_colors,
    ifelse(is.list(df), NROW(df[[2]]), NROW(df)) / length(var_colors) + 1
  )
  var_shapes <- rep(
    var_shapes,
    ifelse(is.list(df), NROW(df[[2]]), NROW(df)) / length(var_shapes) + 1
  )
  AVE_colors <- rep(AVE_colors, NROW(df) / length(AVE_colors) + 1)

  # Call plotting function
  p <- plot_function(
    df = df, title = title, x = x, block = block, comp = comp,
    theme_RGCCA = theme_RGCCA, cex_main = cex_main, cex_sub = cex_sub,
    cex_point = cex_point, sample_colors = sample_colors,
    sample_shapes = sample_shapes, show_labels = show_labels,
    repel = repel, var_colors = var_colors, var_shapes = var_shapes,
    AVE_colors = AVE_colors, show_arrows = show_arrows
  )
  if (is(p, "ggplot")) plot(p, ...)
  invisible(p)
}
