#' Check blocks
#'
#' check_blocks runs several checks on the blocks and transform them in
#' order to ensure that the blocks can be analysed properly.
#'
#' check_blocks performs the following checks and apply the following
#' transformations to the blocks:
#' \itemize{
#'   \item If a single block is given as a data frame or an array, \code{blocks}
#'   is transformed into a list with the block as its unique element. Otherwise,
#'   if \code{blocks} is not a list, an error is raised.
#'   \item Coerce each element of \code{blocks} to an array.
#'   \item Make sure that all the blocks apart from the response block are
#'   quantitative.
#'   \item Add missing names to \code{blocks}.
#'   \item Add missing dimnames to each block and prefix dimnames with
#'   block names if some dimnames are duplicated between blocks.
#'   \item Check blocks' primary dimnames. Raises an error if a block has
#'   duplicated primary dimnames. Several scenario are possible:
#'   \itemize{
#'     \item If all blocks are missing primary dimnames,  primary dimnames are
#'     created.
#'     \item If a block is missing  primary dimnames and all other blocks'
#'     primary dimnames match, missing  primary dimnames are copied from the
#'     other blocks.
#'     \item If a block is missing  primary dimnames but other blocks' have none
#'     matching  primary dimnames, an error is raised.
#'   }
#'   \item If blocks have different number of variables on the primary
#'   dimension, fibers filled with NA values are added to the blocks with
#'   missing variables on the primary dimension. Blocks' variables on the
#'   primary dimension are permuted so that every block has the same primary
#'   dimnames in the same order.
#' }
#' @inheritParams rgcca
#' @importFrom abind abind
#' @noRd
check_blocks <- function(blocks, quiet = FALSE, response = NULL) {
  blocks <- check_blocks_is_list(blocks)
  blocks <- check_blocks_data_structure(blocks)
  blocks <- check_blocks_quantitative(blocks, response)
  blocks <- check_blocks_names(blocks)
  blocks <- check_blocks_dimnames(blocks, 1, quiet)
  blocks <- check_blocks_align(blocks, 1)

  invisible(blocks)
}

check_blocks_is_list <- function(blocks) {
  # Check that there is either a single block or a list of blocks
  if (is.array(blocks) || is.data.frame(blocks)) blocks <- list(blocks)
  if (!is.list(blocks)) stop_rgcca(paste("blocks must be a list."))
  return(blocks)
}

check_blocks_data_structure <- function(blocks) {
  blocks <- lapply(blocks, function(x) {
    if (is.array(x)) {
      return(x)
    }
    if (is.data.frame(x)) {
      return(data.matrix(x))
    }
    names_x <- names(x)
    x <- data.matrix(x)
    rownames(x) <- names_x
    return(x)
  })
  return(blocks)
}

check_blocks_quantitative <- function(blocks, response = NULL) {
  response <- ifelse(is.null(response), length(blocks) + 1, response)
  lapply(seq_along(blocks), function(j) {
    x <- blocks[[j]]
    qualitative <- is.character(x) || is.factor(x)
    if (j == response) {
      if (qualitative && (NCOL(x) > 1)) {
        stop_rgcca(
          "unsupported multivariate qualitative block. Block ", j,
          " is a multivariate qualitative block. The method ",
          "is not able to cope with it."
        )
      }
    } else {
      if (qualitative) {
        stop_rgcca(
          "unsupported qualitative block. Block ", j,
          " is a qualitative block but is not the response block. The method ",
          "is not able to cope with it."
        )
      }
    }
  })
  return(blocks)
}

check_blocks_names <- function(blocks) {
  # Add block names if some are missing
  if (is.null(names(blocks))) names(blocks) <- rep("", length(blocks))
  for (j in which(names(blocks) == "")) {
    names(blocks)[j] <- paste0("block", j)
  }
  return(blocks)
}

check_blocks_dimnames <- function(blocks, primary = 1, quiet = FALSE) {
  # Create dimnames if missing
  missing_dimnames <- vapply(blocks, function(x) {
    is.null(dimnames(x))
  }, FUN.VALUE = logical(1L))
  blocks[missing_dimnames] <- lapply(blocks[missing_dimnames], function(x) {
    dimnames(x) <- list(NULL)
    return(x)
  })

  # Check dimnames on primary dimension
  blocks <- check_blocks_primary_dimnames(blocks, m = primary)

  # Check dimnames on other dimensions
  blocks <- Map(function(block, name) {
    check_block_secondary_dimnames(block, name, primary, quiet)
  }, blocks, names(blocks))

  # Check for duplicated dimnames across blocks (except primary dim)
  if (any(duplicated(unlist(
    lapply(blocks, function(x) dimnames(x)[-primary])
  )))) {
    if (!quiet) message(
      "Duplicated dimnames across blocks are modified to avoid confusion. \n"
    )

    # Add block name as a prefix to avoid confusion
    blocks[seq_along(blocks)] <- lapply(seq_along(blocks), function(j) {
      block <- blocks[[j]]
      dimnames(block)[-primary] <- lapply(
        seq_along(dim(block))[-primary], function(m) {
          dimnames(block)[[m]] <- paste(
            names(blocks)[j], dimnames(blocks[[j]])[[m]], sep = "_"
          )
        }
      )
      return(block)
    })
  }
  return(blocks)
}

check_blocks_primary_dimnames <- function(blocks, m = 1) {
  # Raise error if dimension m does not exist in one of the blocks
  lapply(seq_along(blocks), function(j) {
    if (m > length(dim(blocks[[j]]))) stop_rgcca(
      "wrong number of dimensions. Dimension ", m, " is the dimension shared ",
      "by all the blocks but is missing for block ", j, "."
    )
  })

  # Raise error if there are duplicated dimnames
  lapply(blocks, function(x) {
    duplicated_names <-
      !is.null(dimnames(x)[[m]]) && any(duplicated(dimnames(x)[[m]]))
    if (duplicated_names) stop_rgcca(
      "blocks have duplicated names on dimension ", m, "."
    )
  })

  empty_names <- vapply(
    blocks, function(x) is.null(dimnames(x)[[m]]), FUN.VALUE = logical(1L)
  )

  # Create dimnames for all blocks if all missing
  if (all(empty_names)) {
    blocks <- lapply(blocks, function(x) {
      dimnames(x)[[m]] <- paste0("S", seq_len(dim(x)[[m]]))
      return(x)
    })
  } else if (any(empty_names)) {
    # If at least one block does not have dimnames, 2 cases arise:
    #   - if all blocks with names have the same dimnames, in the same order,
    #     we fill the missing dimnames with the dimnames of the other blocks;
    #   - otherwise we raise an error.
    matrix_of_names <- do.call(
      cbind, lapply(blocks[!empty_names], function(x) dimnames(x)[[m]])
    )
    is_valid <- all(
      apply(matrix_of_names, 2, identical, matrix_of_names[, 1])
    )
    if (is_valid) {
      blocks[empty_names] <- lapply(blocks[empty_names], function(x) {
        dimnames(x)[[m]] <- matrix_of_names[, 1]
        return(x)
      })
    } else {
      stop_rgcca(
        "some blocks are missing names on dimension ", m, ", and the other ",
        "blocks' names on dimension ", m, " are not consistent."
      )
    }
  }

  return(blocks)
}

check_block_secondary_dimnames <- function(x, n, primary = 1, quiet = FALSE) {
  n_dim <- length(dim(x)[-primary])
  # Set default dimnames if missing
  dimnames(x)[-primary] <- lapply(seq_along(dim(x))[-primary], function(m) {
    if (is.null(dimnames(x)[[m]])) {
      if (n_dim == 1) {
        if (dim(x)[m] == 1) {
          return(n)
        } else {
          return(paste(n, seq_len(dim(x)[m]), sep = "_"))
        }
      } else {
        return(paste(n, m, seq_len(dim(x)[m]), sep = "_"))
      }
    } else {
      return(dimnames(x)[[m]])
    }
  })

  # Add mode and variable number as a suffix if there are duplicates
  duplicated_names <- duplicated(unlist(dimnames(x)[-primary]))
  if (any(duplicated_names)) {
    dimnames(x)[-primary] <- lapply(seq_along(dim(x))[-primary], function(m) {
      if (n_dim == 1) {
        return(paste(dimnames(x)[[m]], seq_len(dim(x)[m]), sep = "_"))
      } else {
        return(paste(dimnames(x)[[m]], m, seq_len(dim(x)[m]), sep = "_"))
      }
    })
    if (!quiet) message(
      "Duplicated dimnames are modified to avoid confusion. \n"
    )
  }

  return(x)
}

check_blocks_align <- function(blocks, m = 1) {
  # Construct union of names on dimension m
  all_names <- Reduce(union, lapply(blocks, function(x) dimnames(x)[[m]]))

  # If one block doesn't have as many fibers on mode m as there
  # are names in all_names, we complete the blocks by
  # adding fibers full of NA on mode m.
  lacking_values <- vapply(blocks, function(x) {
    dim(x)[[m]] != length(all_names)
  }, FUN.VALUE = logical(1L))

  blocks[lacking_values] <- lapply(blocks[lacking_values], function(x) {
    missing_names <- setdiff(all_names, dimnames(x)[[m]])
    extra_dim <- dim(x)
    extra_dim[m] <- length(missing_names)
    y <- array(NA, dim = extra_dim)
    dimnames(y)[[m]] <- missing_names
    rownames(y) <- missing_names
    return(abind(x, y, along = m))
  })

  # Align blocks using names on dimension m
  blocks <- lapply(blocks, function(x) {
    perm <- seq_along(dim(x))
    perm[c(1, m)] <- c(m, 1)
    x <- aperm(x, perm)
    x <- subset_block_rows(x, all_names, drop = FALSE)
    aperm(x, perm)
  })
  return(blocks)
}
