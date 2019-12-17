warn_connection <- function(x)
        warning(
            paste0(
                "By using a ",
                x,
                ", all blocks are connected to this block in the connection matrix and the connection file is ignored."
            )
        )
