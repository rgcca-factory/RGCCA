check_blockx <- function(x, y, blocks){
    x <- check_min_integer(x, y, blocks)
    if (y > length(blocks))
        stop_rgcca(
            paste0(
                x,
                " should be lower than ",
                length(blocks),
                " (the maximum number of blocks), not be equal to ",
                y,
                "."
            ),
            exit_code = 133
        )
    return(x)
}
