# y <- runif(6)
# check_integer("y", y, "vector", T, min = 0)
# y <- matrix(runif(6, 1, 10), 2, 3)
# check_integer("y", y, "matrix", T)
# check_integer("y", y, "matrix")
# check_integer(NA)
# check_integer(.1)
# check_integer(c(1:2))
# check_integer("x", c(0, 0), type = "vector")
check_integer <- function(x, y = x, type = "scalar", float = FALSE, min = 1) {

    if (is.null(y))
        y <- x
    
    if (type %in% c("matrix", "data.frame"))
        y_temp <- y

    y <- tryCatch(
        as.double(as.matrix(y)),
        warning = function(w)
            stop(paste(x, "should be numeric."))
        )

    if (any(is.na(y)))
        stop(paste(x, "should not be NA."))

    if (!is(y, "numeric"))
        stop(paste(x, "should be numeric."))
    
    if (type == "scalar" && length(y) != 1)
        stop(paste(x, "should be of length 1."))

    if (!float)
        y <- as.integer(y)
    
    if (all(y < min))
        stop(paste0(x, " should be higher than or equal to ", min, "."))

    if (type %in% c("matrix", "data.frame"))
        y <- matrix(
            y, 
            dim(y_temp)[1], 
            dim(y_temp)[2],
            dimnames = dimnames(y_temp)
        )

    if (type == "data.frame")
        as.data.frame(y)

    return(y)
}
