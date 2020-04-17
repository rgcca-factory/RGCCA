check_ncomp <- function(ncomp, blocks, min = 1) {
    
    ncomp <- elongate_arg(ncomp, blocks)
    ncomp <- sapply(
        seq(length(ncomp)),
        function(x){
            y <- check_integer("ncomp", ncomp[x], min = min)
            if (y > NCOL(blocks[[x]])) {
                stop(
                    paste0(
                        "ncomp should be comprise between ", min ," and ",
                        NCOL(blocks[[x]]),
                        ", the number of variables of the block (currently equals to ",
                        y,
                        ")."
                    ),
                    exit_code = 126
                )
            }
            else
                return(y)
        }
    )
    check_size_blocks(blocks, "ncomp", ncomp)
    return(ncomp)
}
