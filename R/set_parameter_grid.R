#' Set parameter grid
#'
#' Produce a grid of parameters for rgcca (tau, sparsity or ncomp) that will
#' be evaluated either using cross validation or permutation.
#' @inheritParams rgcca_cv
#' @noRd
set_parameter_grid <- function(par_type, par_length, par_value, blocks,
                               response = NULL, superblock = FALSE) {
  ### Auxiliary functions
  check_param_type <- function(par_value, blocks) {
    is_valid_type <- is.null(par_value) || is.vector(par_value) ||
      (length(dim(par_value)) == 2)
    if (!is_valid_type) {
      stop_rgcca(
        "wrong type of input. par_value must be one of the ",
        "following: NULL, a vector, a matrix or a dataframe."
      )
    }
    is_valid_shape <- (NCOL(par_value) == 1) ||
      (NCOL(par_value) == J) ||
      ((NCOL(par_value) == J + 1) && superblock)
    if (!is_valid_shape) {
      stop_rgcca(
        "wrong shape. If par_value is a matrix or a dataframe,",
        "it must have as many columns as there are blocks (i.e. ",
        length(blocks), ")."
      )
    }
  }

  set_response_value <- function(par_value, response_value) {
    if (is.null(response_value)) {
      return(par_value)
    }
    par_value <- t(apply(par_value, 1, function(x) {
      x[response] <- response_value(x)
      return(x)
    }))
  }

  set_grid <- function(check_function, min_values, max_values,
                       response_value = NULL) {
    # If par_value is null, we generate a matrix with par_length rows
    # by taking values uniformly spaced between the min of possible
    # values and the max of possible values for each block.
    if (is.null(par_value)) {
      par_value <- lapply(seq_along(blocks), function(j) {
        seq(max_values, min_values[j], length.out = par_length)
      })
      par_value <- do.call(cbind, par_value)
      par_value <- set_response_value(par_value, response_value)
      return(list(par_type = par_type, par_value = par_value))
    }
    # If par_value is a vector, we aim to create a matrix out of this
    # vector. Hence we have to check beforehand that par_value is a vector
    # of valid numbers.
    if (is.vector(par_value)) {
      par_value <- check_function(par_value)
      par_value <- lapply(seq_along(par_value), function(j) {
        seq(par_value[j], min_values[j], length.out = par_length)
      })
      par_value <- do.call(cbind, par_value)
      par_value <- set_response_value(par_value, response_value)
      return(list(par_type = par_type, par_value = par_value))
    }
    # If par_value is already a grid, we just check that it is valid.
    par_value <- t(vapply(seq_len(NROW(par_value)), function(i) {
      check_function(par_value[i, ])
    }, FUN.VALUE = double(ncol(par_value))))
    par_value <- set_response_value(par_value, response_value)
    return(list(par_type = par_type, par_value = par_value))
  }

  ### Main function
  J <- length(blocks)
  check_param_type(par_value, blocks)
  ncols <- vapply(blocks, NCOL, FUN.VALUE = integer(1))

  switch(par_type,
    "ncomp" = {
      if (!is.null(response)) ncols <- ncols[-response]
      min_values <- rep(1, J + 1)
      max_values <- min(
        ifelse(superblock, sum(ncols), min(ncols)), par_length
      )
      response_value <- function(x) {
        return(max(x[-response]))
      }
      check_function <- function(x) {
        check_ncomp(x, blocks, response = response, superblock = superblock)
      }
    },
    "tau" = {
      min_values <- rep(0, J + 1)
      max_values <- 1
      response_value <- NULL
      check_function <- function(x) {
        check_penalty(x, blocks, method = "rgcca", superblock = superblock)
      }
    },
    "sparsity" = {
      min_values <- c(1 / sqrt(ncols), 1 / sqrt(sum(ncols)))
      max_values <- 1
      response_value <- NULL
      check_function <- function(x) {
        check_penalty(x, blocks, method = "sgcca", superblock = superblock)
      }
    }
  )
  if (is.null(response)) response_value <- NULL

  param <- set_grid(check_function, min_values, max_values, response_value)

  if (par_type == "ncomp") param$par_value <- round(param$par_value)
  param$par_value <-
    param$par_value[!duplicated(param$par_value), , drop = FALSE]

  return(param)
}
