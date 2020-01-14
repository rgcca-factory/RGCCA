# check_boolean ("x", c(T, F), "vector")
# check_boolean ("x")
# check_boolean ("x", c(T, F))
# check_boolean (NA)
check_boolean <- function(x, y = x, type = "scalar") {

    if (is.null(y))
        y <- x

    if (type == "scalar")
        x = ""

    if (any(is.na(y)))
        stop(paste(x, "should not be NA."))

    if (!is(y, "logical"))
        stop(paste(x, "should be TRUE or FALSE"))

    if (type == "scalar" && length(y) != 1)
        stop(paste(x, "should be of length 1."))
}
