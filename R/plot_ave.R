#' Histogram of Average Variance Explained
#'
#' Histogram of the model quality (based on Average Variance Explained)
#' for each blocks and sorted in decreasing order
#' @inheritParams plot_ind
#' @inheritParams plot_histogram
#' @seealso \code{\link[RGCCA]{rgccad}}, \code{\link[RGCCA]{sgcca}}
#' @export
#' @importFrom ggplot2 ggplot
plot_ave <- function(rgcca_res,
                     cex = 1,
                     title = "Average Variance Explained",
                     colors = NULL,
                     ...) {
  stopifnot(is(rgcca_res, "rgcca"))
  check_integer("cex", cex)

  if (tolower(rgcca_res$call$method) == "pca") {
    rgcca_res$AVE$AVE_X <- rgcca_res$AVE$AVE_X[1]
    rgcca_res$call$ncomp <- rgcca_res$call$ncomp[1]
    rgcca_res$a <- rgcca_res$a[1]
  }

  names(rgcca_res$AVE$AVE_X) <- NULL
  ave <- 100 * unlist(rgcca_res$AVE$AVE_X)
  blocks <- factor(unlist(lapply(
    seq(length(names(rgcca_res$a))),
    function(x) rep(names(rgcca_res$a)[x], rgcca_res$call$ncomp[x])
  )),
  levels = names(rgcca_res$a)
  )

  if (is.null(names(ave))) {
    names(ave) <- rep(1, length(ave))
  }
  ncomp <- as.factor(names(ave))

  y_ave_cum <- lapply(
    lapply(
      rgcca_res$AVE$AVE_X,
      function(x) round(100 * cumsum(x), 1)
    ),
    function(x) c(0, x)
  )
  y_ave_cum <- unlist(lapply(y_ave_cum, function(x) {
    unlist(lapply(
      seq(length(x)),
      function(i) (x[i - 1] + x[i]) / 2
    ))
  }))

  ave_label <- unlist(lapply(rgcca_res$AVE$AVE_X, function(x) {
    round(100 * x, 1)
  }))
  ave_label[ave_label < max(y_ave_cum) / 20] <- ""

  df <- data.frame(ave, blocks, ncomp, stringsAsFactors = FALSE)
  class(df) <- c(class(df), "d_ave")

  p <- ggplot(
    data = df,
    aes(
      x = blocks,
      y = ave,
      fill = ncomp,
      label = ave_label
    )
  )

  p <- plot_histogram(
    p,
    df,
    title,
    cex = cex,
    ...
  ) +
    scale_fill_manual(
      values = color_group(levels(df$ncomp), colors = colors),
      labels = gsub("comp", " ", levels(df$ncomp))
    ) +
    geom_col(position = position_stack(reverse = TRUE)) +
    labs(subtitle = print_comp(rgcca_res, outer = TRUE)) +
    geom_text(aes(y = y_ave_cum), cex = 3.5 * cex, color = "white") +
    labs(fill = "Components")

  return(p)
}
