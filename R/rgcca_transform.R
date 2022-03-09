#' RGCCA Transform
#'
#' Project blocks of data on canonical components extracted with a RGCCA model.
#'
#' @param rgcca_res A fitted RGCCA object (see  \code{\link[RGCCA]{rgcca}}).
#' @param X A list of either dataframes or matrices to be projected.
#' @param X_scaled A boolean indicating if the blocks in X have already been
#' scaled.
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
#' projection <- rgcca_transform(fit.rgcca, X, X_scaled = FALSE)
#' @export
rgcca_transform <- function(rgcca_res, X, X_scaled = TRUE) {
  ### Auxiliary function
  scl_fun <- function(data, center, scale) {
    # Use the scaling parameter of the training set on the new set
    if (length(center) != 0) {
      if (is.null(dim(data))) {
        # Case of data is a vector
        data <- as.matrix(data)
      }
      if (is.null(scale)) scale <- FALSE
      res <- scale(data, center, scale)
    } else {
      res <- data
    }
    return(res)
  }

  ### Check input parameters
  stopifnot(is(rgcca_res, "rgcca"))
  check_boolean("X_scaled", X_scaled)
  if (is.null(names(X))) stop_rgcca("Please provide names for the blocks X.")

  ### Align training blocks and X
  if (!all(names(X) %in% names(rgcca_res$call$blocks))) {
    stop_rgcca(paste0(
      "At least one block from X was not found in the training",
      " blocks. Please check block names."
    ))
  }
  X_train <- rgcca_res$call$blocks[names(X)]
  X <- lapply(seq_along(X), function(j) {
    x <- X[[j]]
    if ((is.null(dim(x)) && !is.null(dim(X_train[[j]]))) ||
      (!is.null(dim(x)) && is.null(dim(X_train[[j]]))) ||
      any(dim(x)[-1] != dim(X_train[[j]])[-1])
    ) {
      stop_rgcca(paste0(
        "Dimensions of blocks do not match for block",
        names(X)[[j]]
      ))
    }
    if (!is.null(dim(x))) x <- x[, colnames(X_train[[j]]), drop = FALSE]
    return(x)
  })

  ### Scale X if needed
  if (!X_scaled) {
    X <- lapply(seq_along(X), function(j) {
      scl_fun(
        X[[j]],
        attr(X_train[[j]], "scaled:center"),
        attr(X_train[[j]], "scaled:scale")
      )
    })
  }

  ### Project X on the space computed using RGCCA
  astar <- rgcca_res$astar[names(X_train)]
  projection <- lapply(seq_along(X), function(j) {
    x <- pm(as.matrix(X[[j]]), astar[[j]])
    rownames(x) <- rownames(X[[j]])
    colnames(x) <- colnames(astar[[j]])
    return(x)
  })
  names(projection) <- names(X_train)

  return(projection)
}
