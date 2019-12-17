check_lower_blocks <- function(x, y, blocks)
    if (y > length(blocks))
        stop(
            paste0(
                x,
                " must be lower than ",
                length(blocks),
                " (the maximum number of blocks), not be equal to ",
                y,
                "."
            ),
            exit_code = 133
        )
