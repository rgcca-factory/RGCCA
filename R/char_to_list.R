#' Convert a character in a vector
#'
#' @param x A character separated by comma
#' @return A vector of characters whitout spaces
#' @examples
#' s = '1,2, 3'
#' char_to_list(s)
#' @export
char_to_list <- function(x) strsplit(gsub(" ", "", as.character(x)), ",")[[1]]
