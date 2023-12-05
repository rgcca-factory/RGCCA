#' Check blocks
#'
#' check_blocks runs several checks on the blocks and transform them in
#' order to ensure that the blocks can be analysed properly.
#'
#' check_blocks performs the following checks and apply the following
#' transformations to the blocks:
#' \itemize{
#'   \item If a single block is given as a data frame or a matrix, \code{blocks}
#'   is transformed into a list with the block as its unique element. Otherwise,
#'   if \code{blocks} is not a list, an error is raised.
#'   \item Coerce each element of \code{blocks} to a matrix.
#'   \item Make sure that all the blocks apart from the response block are
#'   quantitative.
#'   \item Add missing names to \code{blocks}.
#'   \item Add missing column names to each block and prefix column names with
#'   block names if some column names are duplicated between blocks.
#'   \item Check blocks' row names. Raises an error if a block has duplicated
#'   row names. Several scenario are possible:
#'   \itemize{
#'     \item If all blocks are missing row names, row names are created if
#'     \code{allow_unnames} is TRUE, otherwise an error is raised.
#'     \item If a block is missing row names and all other blocks' row names
#'     match, missing row names are copied from the other blocks.
#'     \item If a block is missing row names but other blocks' have none
#'     matching row names, an error is raised.
#'   }
#'   \item If \code{add_NAlines} is FALSE and blocks have different number of
#'   rows, an error is raised. Otherwise, lines filled with NA values are added
#'   to the blocks with missing rows. Blocks' rows are permuted so that every
#'   block has the same row names in the same order.
#' }
#' @inheritParams rgcca
#' @param add_Nalines logical, if TRUE, lines filled with NA are added to blocks
#' with missing rows
#' @param allow_unnames logical, if FALSE, an error is raised if blocks do not
#' have row names
#' @importFrom stats setNames
#' @noRd
check_blocks_mg <- function(blocks, add_NAlines = FALSE, allow_unnames = TRUE, #EG
                         quiet = FALSE, response = NULL, groups = NULL) { #EG
  blocks <- check_blocks_is_list(blocks)
  blocks <- check_blocks_matrix(blocks)
  blocks <- check_blocks_quantitative(blocks, response)
  blocks <- check_blocks_names(blocks, quiet)
  blocks <- check_blocks_colnames(blocks, quiet)
  blocks <- check_blocks_rownames(blocks, allow_unnames, quiet)
  blocks <- check_blocks_align(blocks, add_NAlines, quiet)
  blocks <- check_blocks_groups_mg(blocks, groups) #EG
  
  invisible(blocks)
}

check_blocks_is_list <- function(blocks) {
  # Check that there is either a single block or a list of blocks
  if (is.matrix(blocks) || is.data.frame(blocks)) blocks <- list(blocks)
  if (!is.list(blocks)) stop_rgcca(paste("blocks must be a list."))
  return(blocks)
}

check_blocks_matrix <- function(blocks) {
  blocks <- lapply(blocks, function(x) {
    if (is.matrix(x)) {
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

check_blocks_names <- function(blocks, quiet = FALSE) {
  # Add block names if some are missing
  renamed <- FALSE
  if (is.null(names(blocks))) names(blocks) <- rep("", length(blocks))
  for (x in which(names(blocks) == "")) {
    names(blocks)[x] <- paste0("block", x)
    renamed <- TRUE
  }
  if (!quiet && renamed) {
    message("Missing block names are automatically labeled.")
  }
  return(blocks)
}

check_blocks_colnames <- function(blocks, quiet = FALSE) {
  # Check for empty colnames
  if (any(vapply(
    blocks, function(x) is.null(colnames(x)),
    FUN.VALUE = logical(1)
  ))) {
    if (!quiet) message("Missing colnames are automatically labeled.")
    blocks <- lapply(
      setNames(seq_along(blocks), names(blocks)),
      function(x) {
        block <- blocks[[x]]
        if (is.null(colnames(block))) {
          if (NCOL(block) == 1) {
            colnames(block) <- names(blocks)[x]
          } else {
            colnames(block) <- paste0("V", x, "_", seq_len(NCOL(block)))
          }
        }
        return(block)
      }
    )
  }
  
  # Check for duplicated colnames
  if (any(duplicated(unlist(lapply(blocks, colnames))))) {
    if (!quiet) message("Duplicated colnames are modified to avoid confusion.")
    
    blocks <- lapply(
      setNames(seq_along(blocks), names(blocks)),
      function(x) {
        block <- blocks[[x]]
        colnames(block) <- paste(names(blocks)[x],
                                 colnames(blocks[[x]]),
                                 sep = "_"
        )
        return(block)
      }
    )
  }
  return(blocks)
}

check_blocks_rownames <- function(blocks, allow_unnames = TRUE, quiet = FALSE) {
  # Raise error if duplicated rownames
  lapply(blocks, function(x) {
    if (!is.null(row.names(x)) && any(duplicated(row.names(x)))) {
      stop_rgcca(
        "blocks have duplicated rownames."
      )
    }
  })
  
  # Create rownames for all blocks if all missing
  if (all(vapply(
    blocks, function(x) is.null(row.names(x)),
    FUN.VALUE = logical(1)
  ))) {
    if (allow_unnames) {
      blocks <- lapply(
        blocks,
        function(x) {
          rownames(x) <- paste0("S", seq_len(NROW(x)))
          return(x)
        }
      )
      if (!quiet) message("Missing rownames are automatically labeled.")
    } else {
      stop_rgcca(paste("blocks must have rownames."))
    }
  }
  
  # If at least one block does not have rownames, 2 cases arise:
  #   - if all blocks with names have the same rownames, in the same order,
  #     we fill the missing rownames with the rownames of the other blocks;
  #   - otherwise we raise an error.
  if (any(
    vapply(blocks, function(x) is.null(row.names(x)), FUN.VALUE = logical(1))
  )) {
    matrix_of_rownames <- Reduce(cbind, lapply(blocks, row.names))
    is_valid <- all(
      apply(matrix_of_rownames, 2, identical, matrix_of_rownames[, 1])
    )
    if (is_valid) {
      blocks <- lapply(
        blocks,
        function(x) {
          row.names(x) <- matrix_of_rownames[, 1]
          return(x)
        }
      )
      if (!quiet) message("Missing rownames are automatically labeled.")
    } else {
      stop_rgcca(
        "some blocks are missing rownames, and the other blocks' ",
        "rownames are not consistent."
      )
    }
  }
  return(blocks)
}

check_blocks_align <- function(blocks, add_NAlines = FALSE, quiet = FALSE) {
  # Construct union of rownames
  all_names <- Reduce(union, lapply(blocks, row.names))
  
  # If add_NAlines is FALSE and one block doesn't have as many rows as there
  # are names in all_names, we stop. Otherwise we complete the blocks by
  # adding rows full of NA.
  
  if (any(vapply(blocks, nrow, FUN.VALUE = integer(1)) != length(all_names))) {
    if (add_NAlines) {
      blocks <- lapply(blocks, function(x) {
        missing_names <- setdiff(all_names, rownames(x))
        y <- matrix(NA, nrow = length(missing_names), ncol = ncol(x))
        rownames(y) <- missing_names
        return(rbind(x, y))
      })
    } else {
      stop_rgcca("blocks must have the same rownames.")
    }
  }
  
  # Align blocks using rownames
  blocks <- lapply(
    blocks, function(x) x[row.names(blocks[[1]]), , drop = FALSE]
  )
  return(blocks)
}

check_blocks_groups_mg <- function(blocks, groups = NULL) { #EG whole function
  if (!is.null(groups)) {
    # Check that the multigroup parameter groups is a factor of the right length
    groups <- droplevels(as.factor(groups))
    if (length(groups) != nrow(blocks[[1]])) {
      stop_rgcca(paste("groups should be a factor of length ", NROW(blocks[[1]]),"."))
    }
    
    # Check that each group has at least 2 representatives
    if (any(table(groups) < 2)) {
      stop_rgcca(paste("All groups should have more than one representative."))
    }
    
    # Check that only one block was provided
    if (length(blocks) > 1) {
      stop_rgcca(paste("Please provide only one block to use the multi-group framework."))
    }
    
    # Check for NAs in groups
    if (any(is.na(groups))) {
      stop_rgcca(paste("Argument groups cannot accept NAs."))
    }
    
    # Check for NA rows in blocks
    if (any(sapply(blocks, function(X) {
      any(apply(X, 1, function(x) all(is.na(x))))}))) {
      # TRUE if at least one block has at least one row with all missing values
      missing_rows_list <- lapply(blocks, function(X) {
        apply(X, 1, function(x) {all(is.na(x))})
      })
      blocks[seq_along(blocks)] <- lapply(seq_along(blocks), function(j) {
        blocks[[j]][!missing_rows_list[[j]],]})
      if (!quiet) {
        message(paste('Rows that had all values missing were removed.'))
      }
    }
    
    # Check for NA columns in blocks
    if (any(sapply(blocks, function(X) {
      any(apply(X, 2, function(x) all(is.na(x))))}))) {
      # TRUE if at least one block has at least one column with all missing values
      missing_cols_list <- lapply(blocks, function(X) {
        apply(X, 2, function(x) {all(is.na(x))})
      })
      blocks[seq_along(blocks)] <- lapply(seq_along(blocks), function(j) {
        blocks[[j]][,!missing_cols_list[[j]]]})
      if (!quiet) {
        message(paste('Columns that had all values missing were removed.'))
      }
    }
    
    # Split blocks into groups
    blocks <- split(x = data.frame(blocks[[1]]), # split only works on data.frames and vectors (not matrices)
                    f = groups, drop = TRUE) # drop=T to drop levels that do not occur
    blocks <- lapply(blocks, data.matrix)
  }
}
