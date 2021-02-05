scaling <- function(blocks, scale = TRUE, bias = TRUE, scale_block = TRUE) {
    if(scale){
        # Standardization of the variables of each block
        blocks <- lapply(blocks,
                         function(x) scale_array(x, scale = TRUE, bias = bias)
                         )

        # Each block is divided by the square root of its number of variables
        if(scale_block){
           blocks <- lapply(blocks,
                            function(x) {y <- x / sqrt(NCOL(x))
                             attr(y, "scaled:scale") <-
                                 attr(x, "scaled:scale")* sqrt(NCOL(x))
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
