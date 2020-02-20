#' Test for character vector
#'
#' Tests if a dataframe is composed only by qualitative variables
#'
#' @param x A matrix or a vector
#' @param type Character for a function ("any" by default) among "any" or "all"
#' @return A bolean for the presence (FALSE) or the absence (TRUE) of at least
#' one quantitative variable

is.character2 <- function(x, type = "any") {
    # is. character() consider a string with '1.2' as a character, not this function.
    # NA are produced by converting a character into an integer as.vector, avoid
    # factors of character in integer without NA
 
    # NA tolerance :
    
    x <- as.vector(x)
    
    get(type)(
        is.na(
            sapply(
                na.omit(x[x != "NA"]), 
                function(i) {
                tryCatch(
                    as.integer(i),
                    warning = function(w) NA
                )
            })))
}

