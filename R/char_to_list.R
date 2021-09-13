# Convert a character in a vector
#
# @param x A character separated by comma
# @return A vector of characters whitout spaces


char_to_list <- function(x) strsplit(gsub(" ", "", as.character(x)), ",")[[1]]
