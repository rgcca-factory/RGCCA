check_spars <- function(blocks, tau, type = "rgcca") {
    # sparsity : A vector of integer giving the spasity parameter for SGCCA (sparsity)
    # Stop the program if at least one sparsity parameter is not in the required interval

    if (tolower(type) == "sgcca") {
        #the minimum value avalaible
        min_sparsity <- lapply(blocks, function(x) 1 / sqrt(NCOL(x)))

        # Check sparsity varying between 1/sqrt(pj) and 1
        tau <- mapply(
            function(x, y) {
            x <- check_integer("sparsity", x, float = TRUE, min = 0)
            if (x < y | x > 1)
                stop(
                    paste0(
                        "Sparsity parameter is equals to ",
                        x,
                        ". For SGCCA, it must be comprise between 1/sqrt(number_column) (i.e., ",
                        toString(unlist(
                            lapply(min_sparsity, function(x)
                                ceiling(x * 100) / 100)
                        ))
                        ,
                        ") and 1."
                    ),
                    exit_code = 132
                )
                else
                    x
        }, tau, min_sparsity)
    }
    
    invisible(tau)
}
