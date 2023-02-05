#' Center and scale a list of blocks
#' @inheritParams rgcca
#' @param na.rm A logical, if TRUE, NA values are replaced by 0 to
#' compute scaling parameters.
#' @noRd
scaling <- function(blocks, scale = TRUE, bias = TRUE,
                    scale_block = "inertia", na.rm = TRUE) {
  if (isTRUE(scale_block)) scale_block <- "inertia"
  sqrt_N <- sqrt(NROW(blocks[[1]]) + bias - 1)

  # Center and eventually scale the blocks
  blocks <- lapply(
    blocks,
    function(x) scale2(x, scale = scale, bias = bias)
  )

  # Scale each block by a constant if requested
  if (scale_block == "lambda1") {
    blocks <- scale_lambda1(blocks, sqrt_N, scale, na.rm = na.rm)
  } else if (scale_block == "inertia") {
    blocks <- scale_inertia(blocks, sqrt_N, scale, na.rm = na.rm)
  }

  return(blocks)
}
