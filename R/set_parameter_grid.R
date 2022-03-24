# Set parameter grid
#
# Produce a grid of parameters for rgcca (tau, sparsity or ncomp) that will
# be evaluated either using cross validation or permutation.
set_parameter_grid <- function(par_type, par_length, par_value, blocks,
                               response) {
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
      (NCOL(par_value) == length(blocks))
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
    return(par_value[!duplicated(par_value), ])
  }

  set_grid <- function(check_function, min_values, max_values, length_values,
                       response_value = NULL) {
    # If par_value is null, we generate a matrix with par_length rows
    # by taking values uniformly spaced between the min of possible
    # values and the max of possible values for each block.
    if (is.null(par_value)) {
      par_value <- lapply(seq_along(blocks), function(j) {
        seq(min_values[j], max_values, length.out = length_values)
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
      par_value <- lapply(seq_along(blocks), function(j) {
        seq(min_values[j], par_value[j], length.out = length_values)
      })
      par_value <- do.call(cbind, par_value)
      par_value <- set_response_value(par_value, response_value)
      return(list(par_type = par_type, par_value = par_value))
    }
    # If par_value is already a grid, we just check that it is valid.
    par_value <- t(vapply(seq(nrow(par_value)), function(i) {
      check_function(par_value[i, ])
    }, FUN.VALUE = double(ncol(par_value))))
    par_value <- set_response_value(par_value, response_value)
    return(list(par_type = par_type, par_value = par_value))
  }

  ### Main function
  check_param_type(par_value, blocks)
  ncols <- vapply(blocks, NCOL, FUN.VALUE = integer(1))

  switch(par_type,
    "ncomp" = {
      min_values <- rep(1, length(blocks))
      max_values <- min(min(ncols[-response]), par_length)
      length_values <- min(min(ncols[-response], par_length))
      response_value <- function(x) {
        return(max(x[-response]))
      }
      check_function <- function(x) {
        check_ncomp(x, blocks, response = response)
      }
    },
    "tau" = {
      min_values <- rep(0, length(blocks))
      max_values <- 1
      length_values <- par_length
      response_value <- function(x) {
        return(x[response])
      }
      check_function <- function(x) {
        check_penalty(x, blocks, method = "rgcca")
      }
    },
    "sparsity" = {
      min_values <- 1 / sqrt(ncols)
      max_values <- 1
      length_values <- par_length
      response_value <- function(x) {
        return(1)
      }
      check_function <- function(x) {
        check_penalty(x, blocks, method = "sgcca")
      }
    }
  )
  set_grid(
    check_function, min_values, max_values, length_values, response_value
  )
}
