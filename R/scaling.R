#' Center and scale a list of blocks
#' @inheritParams rgcca
#' @param na.rm A logical, if TRUE, NA values are replaced by 0 to
#' compute scaling parameters.
#' @noRd
scaling <- function(blocks, scale = TRUE, bias = FALSE,
                    scale_block = "inertia", na.rm = TRUE) {
  if (isTRUE(scale_block)) scale_block <- "inertia"
  sqrt_N <- sqrt(NROW(blocks[[1]]) + bias - 1)

  blocks <- lapply(blocks, function(x) {
    # Store dim and dimnames
    dim_x <- dim(x)
    dimnames_x <- dimnames(x)

    num_dims <- length(dim_x)
    n_var <- dim_x[2]

    ###### CENTRER : by column

    # Unfold the array if needed
    if (length(dim_x) > 2) {
      x <- matrix(x, nrow = nrow(x))
    }

    # Center and eventually scale the blocks
    x <- scale(x, center = TRUE, scale = FALSE)

    # Go back to a tensor
    y <- array(x, dim = dim_x)
    dimnames(y) <- dimnames_x
    ctr <- attr(x, "scaled:center")
    
    ####### SCALE : by slice/variable

    # Unfold the array if needed
    if (num_dims > 2) {
      scale_block = "inertia"

      # Flatten tensor
      mat <- matrix(aperm(y, c(setdiff(seq_len(num_dims), 2), 2)), ncol = n_var)

      # Center and eventually scale the blocks
      mat <- scale_new(mat, scale = scale, bias = bias)
      
      mat <- scale_inertia(mat, sqrt_N, scale, na.rm = na.rm)
      
      # Go back to a tensor
      inv_perm <- match(seq_len(num_dims), c(setdiff(seq_len(num_dims), 2), 2)) #La fonction match(a, b) renvoie, pour chaque élément de a, l’indice où il se trouve dans b
      x_perm <- array(mat, dim = dim_x[c(setdiff(seq_len(num_dims), 2), 2)])
      x_rec <- aperm(x_perm, inv_perm)
      
      y <- x_rec

    } else {
      mat <- x
      mat <- scale_new(mat, scale = scale, bias = bias)
      if (scale_block == "lambda1") {
        mat <- scale_lambda1(mat, sqrt_N, scale, na.rm = na.rm)
      } else if (scale_block == "inertia") {
        mat <- scale_inertia(mat, sqrt_N, scale, na.rm = na.rm)
      }
      y <- mat
    }

    dimnames(y) <- dimnames_x
    
    attr(y, "scaled:center") <- ctr
    attr(y, "scaled:scale") <- attr(mat, "scaled:scale")
    
    return(y)
  })
  
  return(blocks)
}