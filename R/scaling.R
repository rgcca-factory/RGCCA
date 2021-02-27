# TODO: check it works for tensors if scale = FALSE
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
                    attr(y, "scaled:center") <- attr(x, "scaled:center") / sqrt(nb_var)
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
                if(NROW(x) > NCOL(x))
                {
                    covarMat <- cov2(x, bias = bias)
                }
                else
                {
                    covarMat <- 1/N*(x%*%t(x))
                }
                variance_block <- sum(diag(covarMat))
                out <- x / sqrt(variance_block)
                attr(out, "scaled:scale") <-
                    rep(sqrt(sum(diag(covarMat))), NCOL(x))
                return(out)
            })
        }

    }
    return(blocks)
}

# Function to apply scaling/centering to new blocks knowing centering and 
# scaling factors.
apply_scaling <- function(blocks, centering_factors, scaling_factors) {
    lapply(1:length(blocks), function(x) {
        if (length(dim(blocks[[x]])) > 2) {
            block = t(apply(blocks[[x]], 1, function(y) y - centering_factors[[x]]))
            block = apply(block, 2, function(y) y / scaling_factors[[x]])
            return(array(block, dim = dim(blocks[[x]])))
        }
        if (length(dim(blocks)) == 2) {
            block = t(apply(blocks[[x]], 1, function(y) y / scaling_factors[[x]]))
            return(t(apply(block, 1, function(y) y - centering_factors[[x]])))
        }
        return(blocks[[x]] / scaling_factors[[x]] - centering_factors[[x]])
    })
}
