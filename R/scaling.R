scaling <- function(blocks, scale = TRUE, bias = TRUE, scale_block = TRUE) {
    if(scale){
        # Standardization of the variables of each block
        blocks <- lapply(blocks,
                         function(x) scale2(x, scale = TRUE, bias = bias)
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
                             scale2(x, scale = FALSE, bias = bias)
                         )

        if (scale_block) {
            sqrt_N = sqrt(NROW(blocks[[1]]) + bias - 1)
            blocks <- lapply(blocks, function(x) {
                fac <- 1/sqrt_N * norm(x, type = "F")
                out <- x / fac
                attr(out, "scaled:scale") <- rep(fac, NCOL(x))
                return(out)
            })
        }

    }
    return(blocks)
}
