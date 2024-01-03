#' Center and scale a list of blocks
#' @inheritParams rgcca
#' @noRd
scaling <- function(blocks, scale = TRUE, bias = TRUE,
                    scale_block = "inertia", groups = NULL) {
  if (isTRUE(scale_block)) scale_block <- "inertia"
  sqrt_N <- sqrt(NROW(blocks[[1]]) + bias - 1)
  
  if (scale) {
    # Standardization of the variables of each block
    # Or, in the multi-group case, centering and normalization of the variables
    blocks <- lapply(
      blocks,
      function(x) scale2(x, scale = TRUE, bias = bias, groups = groups)
    )
    
    if (!is.null(groups)) return(blocks)
    
    # Each block is divided by a constant that depends on the block
    if (scale_block == "lambda1") {
      blocks <- lapply(blocks, function(x) {
        lambda <- sqrt(ifelse(ncol(x) < nrow(x),
                              eigen(crossprod(x / sqrt_N))$values[1],
                              eigen(tcrossprod(x / sqrt_N))$values[1]
        ))
        y <- x / lambda
        attr(y, "scaled:scale") <- attr(x, "scaled:scale") * lambda
        return(y)
      })
    } else if (scale_block == "inertia") {
      blocks <- lapply(blocks, function(x) {
        y <- x / sqrt(NCOL(x))
        attr(y, "scaled:scale") <- attr(x, "scaled:scale") * sqrt(NCOL(x))
        return(y)
      })
    }
  } else {
    blocks <- lapply(blocks, function(x) {
      scale2(x, scale = FALSE, bias = bias)
    })
    if (scale_block == "lambda1") {
      blocks <- lapply(blocks, function(x) {
        lambda <- sqrt(ifelse(ncol(x) < nrow(x),
                              eigen(crossprod(x / sqrt_N))$values[1],
                              eigen(tcrossprod(x / sqrt_N))$values[1]
        ))
        y <- x / lambda
        attr(y, "scaled:scale") <- rep(lambda, NCOL(x))
        return(y)
      })
    } else if (scale_block == "inertia") {
      blocks <- lapply(blocks, function(x) {
        fac <- 1 / sqrt_N * norm(x, type = "F")
        out <- x / fac
        attr(out, "scaled:scale") <- rep(fac, NCOL(x))
        return(out)
      })
    }
  }
  return(blocks)
}
