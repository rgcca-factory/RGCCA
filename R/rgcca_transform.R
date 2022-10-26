#' RGCCA Transform
#'
#' Project blocks of data on canonical components extracted with a RGCCA model.
#'
#' @param rgcca_res A fitted RGCCA object (see  \code{\link[RGCCA]{rgcca}}).
#' @param blocks_test A list of either dataframes or matrices to be projected.
#' @examples
#' data("Russett")
#' blocks <- list(
#'   agriculture = Russett[, 1:3],
#'   industry = Russett[, 4:5],
#'   politic = Russett[, 6:11]
#' )
#' C <- connection <- 1 - diag(3)
#' A <- lapply(blocks, function(x) x[1:32, ])
#' fit.rgcca <- rgcca(A,
#'   connection = C, tau = c(0.7, 0.8, 0.7),
#'   ncomp = c(3, 2, 4), scale_block = FALSE, superblock = FALSE
#' )
#' X <- lapply(blocks, function(x) x[39:47, ])
#' projection <- rgcca_transform(fit.rgcca, X)
#' @export
rgcca_transform <- function(rgcca_res, blocks_test) {
  ### Auxiliary function
  scl_fun <- function(data, center, scale) {
    # Use the scaling parameter of the training set on the new set
    if (length(center) != 0) {
      if (is.null(scale)) scale <- FALSE
      data <- scale(data, center, scale)
    }
    return(data)
  }

  ### Check input parameters
  stopifnot(is(rgcca_res, "rgcca"))
  if (is.null(names(blocks_test))) {
    stop_rgcca("Please provide names for blocks_test.")
  }

  ### Align training blocks and blocks_test
  if (!all(names(blocks_test) %in% names(rgcca_res$call$blocks))) {
    stop_rgcca(paste0(
      "At least one block from blocks_test was not found in the training",
      " blocks. Please check block names."
    ))
  }
  X_train <- rgcca_res$call$blocks[names(blocks_test)]
  blocks_test <- lapply(seq_along(blocks_test), function(j) {
    x <- as.matrix(blocks_test[[j]])
    y <- as.matrix(X_train[[j]])
    if (!all(colnames(y) %in% colnames(x))) {
      stop_rgcca(
        "Some columns are missing for test block ",
        names(blocks_test)[[j]]
      )
    }
    x <- x[, colnames(y), drop = FALSE]
    return(x)
  })

  ### Scale blocks_test if needed
  blocks_test <- lapply(seq_along(blocks_test), function(j) {
    scl_fun(
      blocks_test[[j]],
      attr(X_train[[j]], "scaled:center"),
      attr(X_train[[j]], "scaled:scale")
    )
  })

  ### Project blocks_test on the space computed using RGCCA
  astar <- rgcca_res$astar[names(X_train)]
  projection <- lapply(seq_along(blocks_test), function(j) {
    x <- pm(as.matrix(blocks_test[[j]]), astar[[j]])
    rownames(x) <- rownames(blocks_test[[j]])
    colnames(x) <- colnames(astar[[j]])
    return(x)
  })
  names(projection) <- names(X_train)

  return(projection)
}
