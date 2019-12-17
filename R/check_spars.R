check_spars <- function(blocks, tau, type = "rgcca") {
    # c1 : A vector of integer giving the spasity parameter for SGCCA (c1)
    # Stop the program if at least one c1 parameter is not in the required interval

    if (tolower(type) == "sgcca") {
        #the minimum value avalaible
        min_c1 <- lapply(blocks, function(x) 1 / sqrt(NCOL(x)))

        # Check c1 varying between 1/sqrt(pj) and 1
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
                            lapply(min_c1, function(x)
                                ceiling(x * 100) / 100)
                        ))
                        ,
                        ") and 1."
                    ),
                    exit_code = 132
                )
                else
                    x
        }, tau, min_c1)
    }
    
    invisible(tau)
}
