elongate_arg <- function(x, blocks) {
    if (length(x) == 1)
        rep(x, length(blocks))
    else
        x
}
