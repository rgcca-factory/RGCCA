#' Reduce dimensionality using RGCCA
#'
#' This function projects testing blocks using the block weight vectors of a
#' fitted RGCCA object.
#'
#' @param rgcca_res A fitted RGCCA object (see  \code{\link[RGCCA]{rgcca}}).
#' @param blocks_test A list of blocks (data.frame or matrix) to be projected.
#' @return A list of matrices containing the projections of the test blocks
#' using the block weight vectors of a fitted RGCCA object.
#' @examples
#' data("Russett")
#' blocks <- list(
#'   agriculture = Russett[, 1:3],
#'   industry = Russett[, 4:5],
#'   politic = Russett[, 6:11])
#'
#' Xtrain <- lapply(blocks, function(x) x[1:32, ])
#' Xtest <- lapply(blocks, function(x) x[33:47, ])
#' fit_rgcca <- rgcca(Xtrain, ncomp = 2)
#' projection <- rgcca_transform(fit_rgcca, Xtest)
#' @export
rgcca_transform <- function(rgcca_res, blocks_test = rgcca_res$call$blocks) {
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
  if (!all(names(blocks_test) %in% names(rgcca_res$blocks))) {
    stop_rgcca(paste0(
      "At least one block from blocks_test was not found in the training",
      " blocks. Please check block names."
    ))
  }
  X_train <- rgcca_res$blocks[names(blocks_test)]
  blocks_test <- lapply(seq_along(blocks_test), function(j) {
    x <- to_mat(blocks_test[[j]])
    y <- to_mat(X_train[[j]])
    # Deal with qualitative block
    if (rgcca_res$opt$disjunction) {
      j_train <- which(names(rgcca_res$blocks) == names(blocks_test)[j])
      if (j_train == rgcca_res$call$response) {
        x <- as_disjunctive(x)
      }
    }
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
  # If there is a superblock with orthogonal components, the superblock
  # is constructed and projected
  if (rgcca_res$call$superblock && rgcca_res$call$comp_orth) {
    superblock_test <- do.call(cbind, blocks_test)
    projection <- list(
      superblock = pm(as.matrix(superblock_test), rgcca_res$astar)
    )
    rownames(projection[[1]]) <- rownames(blocks_test[[1]])
    colnames(projection[[1]]) <- colnames(rgcca_res$astar)
  # Otherwise we directly use astar to project the individual blocks
  } else {
    astar <- rgcca_res$astar[names(X_train)]
    projection <- lapply(seq_along(blocks_test), function(j) {
      x <- pm(as.matrix(blocks_test[[j]]), astar[[j]])
      rownames(x) <- rownames(blocks_test[[j]])
      colnames(x) <- colnames(astar[[j]])
      return(x)
    })
    names(projection) <- names(X_train)
  }
  return(projection)
}
