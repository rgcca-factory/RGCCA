#' Center and scale a list of blocks
#' @inheritParams rgcca
#' @param na.rm A logical, if TRUE, NA values are replaced by 0 to
#' compute scaling parameters.
#' @noRd
scaling <- function(blocks, scale = TRUE, bias = TRUE,
                    scale_block = "inertia", na.rm = TRUE) {
  if (isTRUE(scale_block)) scale_block <- "inertia"
  sqrt_N <- sqrt(NROW(blocks[[1]]) + bias - 1)

  blocks <- lapply(blocks, function(x) {
    # Store dim and dimnames
    dim_x <- dim(x)
    dimnames_x <- dimnames(x)

    # Unfold the array if needed
    if (length(dim_x) > 2) {
      x <- matrix(x, nrow = nrow(x))
    }

    # Center and eventually scale the blocks
    x <- scale2(x, scale = scale, bias = bias)

    # Scale each block by a constant if requested
    if (scale_block == "lambda1") {
      x <- scale_lambda1(x, sqrt_N, scale, na.rm = na.rm)
    } else if (scale_block == "inertia") {
      x <- scale_inertia(x, sqrt_N, scale, na.rm = na.rm)
    }

    # Go back to a tensor
    y <- array(x, dim = dim_x)
    dimnames(y) <- dimnames_x
    attr(y, "scaled:center") <- attr(x, "scaled:center")
    attr(y, "scaled:scale") <- attr(x, "scaled:scale")

    return(y)
  })

  return(blocks)
}
