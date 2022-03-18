#' Plot the two components of a RGCCA
#'
#' Plot the two components of a RGCCA
#'
#' @inheritParams plot2D
#' @param resp A vector of characters corresponding either to a qualitative
#' variable with levels or a continuous variable
#' @param response_name A character giving the legend title
#' @param predicted A list containing as  2nd element a matrix of predicted
#' components
#' @param legend A logical value indicatif if the legend should be plotted
#' @param ... Further graphical parameters (see plot2D functions)
#' @examples
#' coord <- lapply(
#'   seq(3),
#'   function(x) matrix(runif(15 * 2, min = -1), 15, 2)
#' )
#' AVE_X <- lapply(seq(3), function(x) runif(2))
#' for (i in 1:length(coord)) {
#'   row.names(coord[[i]]) <- seq(15)
#' }
#' rgcca_out <- list(
#'   Y = coord, AVE = list(AVE_X = AVE_X),
#'   call = list(blocks = coord, ncomp = rep(2, 3))
#' ) # TODO
#' # Using a superblock
#' resp <- as.matrix(rep(LETTERS[seq(3)], each = 5))
#' row.names(resp) <- seq(15)
#' rgcca_out$call$method <- "rgcca"
#' class(rgcca_out) <- "rgcca"
#' plot_ind(rgcca_out, resp)
#' # Using the first block
#' resp <- as.matrix(runif(15, min = -15, max = 15))
#' row.names(resp) <- seq(15)
#' plot_ind(rgcca_out, resp, 1, 2, 1)
#' data(Russett)
#' X_agric <- as.matrix(Russett[, c("gini", "farm", "rent")])
#' X_ind <- as.matrix(Russett[, c("gnpr", "labo")])
#' X_polit <- as.matrix(Russett[, c("demostab", "dictator")])
#' A <- list(X_agric, X_ind, X_polit)
#' C <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3)
#' result.rgcca <- rgcca(A,
#'   connection = C, tau = c(1, 1, 1), scheme = "factorial",
#'   scale = TRUE, ncomp = rep(2, 3)
#' )
#' plot_ind(result.rgcca, i_block = 1)
#' @export
#' @importFrom ggplot2 ggplot
plot_ind <- function(rgcca_res,
                     resp = rep(1, NROW(rgcca_res$Y[[1]])),
                     compx = 1,
                     compy = 2,
                     i_block = length(rgcca_res$Y),
                     text = TRUE,
                     i_block_y = i_block,
                     response_name = "Response",
                     no_overlap = FALSE,
                     predicted = NULL,
                     title = paste0(
                       names(rgcca_res$call$blocks)[i_block],
                       ": Sample space"
                     ),
                     legend = TRUE,
                     cex = 1,
                     cex_main = 14 * cex,
                     cex_sub = 12 * cex,
                     cex_point = 3 * cex,
                     cex_lab = 10 * cex,
                     ...) {
  if (is.null(i_block_y)) {
    i_block_y <- i_block
  }

  df <- get_comp(
    rgcca_res = rgcca_res,
    resp = resp,
    comps = c(compx, compy),
    i_blocks = c(i_block, i_block_y),
    predicted = predicted
  )
  class(df) <- c(class(df), "d_ind")
  if (!is.null(predicted)) {
    p <- ggplot(df, aes(df[, 1], df[, 2], group = resp, color = resp))
  } else if (length(unique(as.matrix(df$resp))) > 5 &&
    !is.character(as.vector(df$resp))) {
    p <- ggplot(df, aes(df[, 1], df[, 2], group = resp, color = resp))
  } else {
    p <- NULL
  }


  p <- plot2D(
    rgcca_res,
    df,
    title,
    df$resp,
    response_name,
    compx,
    compy,
    i_block,
    p,
    text,
    i_block_y,
    no_overlap = no_overlap,
    cex = cex,
    cex_main = cex_main,
    cex_sub = cex_sub,
    cex_point = cex_point,
    cex_lab = cex_lab,
    ...
  )

  # remove legend if missing
  if (length(unique(df$resp)) == 1 || !legend) {
    p <- p + theme(legend.position = "none")
  }

  return(p)
}
