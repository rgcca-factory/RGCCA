check_min_integer <- function(x , y, blocks, min = 1) {
    x <- check_integer(x, y, min = min)
    check_lower_blocks(x, y, blocks)
    return(x)
}
