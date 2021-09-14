# TODO: check it works for tensors if scale = FALSE
#' @export scaling
#' @export apply_scaling
scaling <- function(blocks, scale = TRUE, bias = TRUE, scale_block = TRUE) {
    if(scale){
        # Standardization of the variables of each block
        blocks <- lapply(blocks,
                         function(x) scale_array(x, scale = TRUE, bias = bias)
                         )

        # Each block is divided by the square root of its number of variables
        if(scale_block){
            blocks <- lapply(blocks, function(x) {
                if (length(dim(x)) > 2) {
                    nb_var <- prod(dim(x)[-1])
                    y <- x / sqrt(nb_var)
                    attr(y, "scaled:scale")  <- attr(x, "scaled:scale") * sqrt(nb_var)
                    attr(y, "scaled:center") <- attr(x, "scaled:center")
                    return(y)
                }
                y <- x / sqrt(NCOL(x))
                attr(y, "scaled:scale") <- attr(x, "scaled:scale")* sqrt(NCOL(x))
                return(y)
        })

        }
    } else{
        blocks <- lapply(blocks,
                         function(x)
                             scale_array(x, scale = FALSE, bias = bias)
                         )

        if (scale_block) {
            N = ifelse(bias, NROW(blocks[[1]]), NROW(blocks[[1]])-1)
            blocks <- lapply(blocks, function(x) {
                DIM = dim(x)
                dn  = dimnames(x)
                if (length(dim(x)) > 2) x = matrix(x, nrow(x))
                s = norm(matrix(x), type = "F") / sqrt(N)
                x = x / s
                if (length(DIM) > 2) {
                    x = array(x, dim = DIM)
                    attr(x, "scaled:scale") <- rep(s, prod(DIM[-1]))
                } else {
                    attr(x, "scaled:scale") <- rep(s, NCOL(x))
                }
                dimnames(x) = dn
                return(x)
            })
        }

    }
    return(blocks)
}

# # Function to apply scaling/centering to new blocks knowing centering and
# # scaling factors.
# apply_scaling <- function(blocks, centering_factors, scaling_factors) {
#     lapply(1:length(blocks), function(x) {
#         if (length(dim(blocks[[x]])) > 2) {
#             block = t(apply(blocks[[x]], 1, function(y) y - centering_factors[[x]]))
#             block = apply(block, 2, function(y) y / scaling_factors[[x]])
#             return(array(block, dim = dim(blocks[[x]])))
#         }
#         if (length(dim(blocks[[x]])) == 2 && dim(blocks[[x]])[2] > 1) {
#             block = t(apply(blocks[[x]], 1, function(y) y - centering_factors[[x]]))
#             return(t(apply(block, 1, function(y) y / scaling_factors[[x]])))
#         }
#         return((blocks[[x]] - centering_factors[[x]]) / scaling_factors[[x]])
#     })
# }

apply_scaling <- function(blocks, centering_factors, scaling_factors) {
    lapply(1:length(blocks), function(x) {
        DIM = dim(blocks[[x]])
        if (length(DIM) > 2) {
            block = t(apply(blocks[[x]], 1, function(y) y - centering_factors[[x]]))
            block = array(block, dim = DIM)
            block = apply(block, -2, function(y) y / scaling_factors[[x]])
            block = aperm(block, c(2, 1, 3:length(DIM)))
            return(block)
        }
        if (length(DIM) == 2 && DIM[2] > 1) {
            block = t(apply(blocks[[x]], 1, function(y) y - centering_factors[[x]]))
            return(t(apply(block, 1, function(y) y / scaling_factors[[x]])))
        }
        return((blocks[[x]] - centering_factors[[x]]) / scaling_factors[[x]])
    })
}
