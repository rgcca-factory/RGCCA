check_compx <- function(x, y, ncomp, blockx) {
    if (y > ncomp[blockx]) {
        stop(
            paste0(
                x,
                " is currently equals to ",
                y,
                " and should be comprise between 1 and ",
                ncomp[blockx],
                " (the number of component for the selected block)."
            ),
            exit_code = 128
        )
    }
    check_integer(x, y, min = 1)
}
