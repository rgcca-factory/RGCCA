#' Test for character vector
#'
#' Tests if a dataframe is composed only by qualitative variables
#'
#' @param x A matrix or a vector
#' @return A bolean for the presence (FALSE) or the absence (TRUE) of at least
#' one quantitative variable
#' @examples
#' x = matrix(c(runif(10), LETTERS[seq(10)]), 10, 2)
#' is.character2(x)
#' # TRUE
#' is.character2(LETTERS[seq(10)])
#' # TRUE
#' @export
is.character2 <- function(x) {
    # is. character() consider a string with '1.2' as a character, not this function.
    # NA are produced by converting a character into an integer as.vector, avoid
    # factors of character in integer without NA
    
    # NA tolerance :
    
    x <- as.vector(x)
    
    any(
        is.na(
            tryCatch(
                as.integer(na.omit(x[x != "NA"])),
                warning = function(w) NA
            )))
}