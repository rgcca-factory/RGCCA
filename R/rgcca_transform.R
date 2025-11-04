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
  ### Auxiliary functions
  scl_fun <- function(data, type, center = NULL, scale = NULL) {
    # Use the scaling parameter of the training set on the new set
    if (type == 'center') {
      if (length(center) != 0) {
        data <- scale(data, center, scale=FALSE)
      }
      return(data)
    }
    
    if (type == 'scale') {
      if (is.null(scale)) scale <- FALSE
      data <- scale(data, center=FALSE, scale)
      return(data)
    }
  }
  
  scl_tens_mat <- function(obj, num_d, n_v, dim_x, dimnames, col_names) {
    # the function enables to matricise a tensor if "obj" is a tensor or retrieve a tensor if "obj" is a matrix
    perm <- c(setdiff(seq_len(num_d), 2), 2)
    
    if (length(dim(obj)) > 2) {
      mat <- matrix(aperm(obj, perm), ncol = n_v)
      colnames(mat) <- col_names
      return(mat)
    } else {
      inv_perm <- match(seq_len(num_d), perm) 
      x_perm <- array(obj, dim = dim_x[perm])
      x_rec <- aperm(x_perm, inv_perm)
      dimnames(x_rec) <- dimnames
      return(x_rec)
    }
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
  names_bl <- names(blocks_test)
  
  ### Center 
  blocks_test <- lapply(seq_along(blocks_test), function(j) {
    # Store dim and dimnames
    dim_x <- dim(blocks_test[[j]])
    dimnames_x <- dimnames(blocks_test[[j]])
    
    # Matricise
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
    
    # Center
    x <- scl_fun(
      x, type = 'center',
      center = attr(X_train[[j]], "scaled:center")
    )
    
    # Go back to a tensor
    x <- array(x, dim = dim_x)
    dimnames(x) <- dimnames_x
    
    return(x)
  })
  
  names(blocks_test) <- names_bl
  
  ### Scale blocks if needed
  blocks_test <- lapply(seq_along(blocks_test), function(j) {
    # Store dim and dimnames
    dim_x <- dim(blocks_test[[j]])
    num_dims <- length(dim_x)
    n_var <- dim_x[2]
    dimnames_x <- dimnames(blocks_test[[j]])
    if (num_dims > 2) {
      x <- scl_tens_mat(blocks_test[[j]], num_dims, n_var, dim_x, 
                        dimnames_x, dimnames_x[[2]])
    } else {
      x <- blocks_test[[j]]
    }
    
    # Scale
    x <- scl_fun(
      x, type = 'scale',
      scale = attr(X_train[[j]], "scaled:scale")
    )
    
    # Go back to tensor if it was a tensor
    if (length(dim_x) > 2) {
      x <- scl_tens_mat(x, num_dims, n_var, dim_x, 
                        dimnames_x, dimnames_x[[2]])
    }
    
    return(x)
  })
  
  names(blocks_test) <- names_bl
  
  ### Matricise
  blocks_test <- lapply(seq_along(blocks_test), function(j) {
    x <- to_mat(blocks_test[[j]])
    return(x)
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
