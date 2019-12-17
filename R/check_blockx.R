check_blockx <- function(x, y, blocks){
    if (y > length(blocks))
        stop(
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
    check_min_integer(x, y, blocks)
}
