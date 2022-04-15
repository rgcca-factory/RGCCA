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
            colnames(block) <- paste0("V", x, "_", seq(NCOL(block)))
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
          rownames(x) <- paste0("S", seq(NROW(x)))
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

check_blocks_character <- function(blocks, no_character = FALSE) {
  # Raise error if characters are present but not allowed
  if (no_character) {
    if (any(vapply(blocks, is.character2, FUN.VALUE = logical(1)))) {
      stop("blocks contain non-numeric values.")
    }

    blocks <- lapply(blocks, as.numeric)
  }
  return(blocks)
}

check_blocks_remove_null_sd <- function(blocks, init = FALSE) {
  if (init) {
    blocks <- remove_null_sd(blocks)$list_m
    for (i in seq_along(blocks)) {
      attributes(blocks[[i]])$nrow <- nrow(blocks[[i]])
    }
  }
  return(blocks)
}

check_blocks <- function(blocks, init = FALSE, n = 2,
                         add_NAlines = FALSE, allow_unnames = TRUE,
                         quiet = FALSE, no_character = FALSE) {
  blocks <- check_blocks_is_list(blocks)
  blocks <- check_blocks_matrix(blocks)
  blocks <- check_blocks_names(blocks, quiet)
  blocks <- check_blocks_colnames(blocks, quiet)
  blocks <- check_blocks_rownames(blocks, allow_unnames, quiet)
  blocks <- check_blocks_remove_null_sd(blocks, init)
  blocks <- check_blocks_align(blocks, add_NAlines, quiet)

  invisible(blocks)
}
