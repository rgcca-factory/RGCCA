#' Keep only the rows with the same names among a list of dataframe
#'
#' @param list_m A list of matrix
#' @return A list of matrix
common_rows <- function(list_m) {
    
    x <- Reduce(intersect, lapply(list_m, row.names))
    
    for (i in seq(length(list_m)))
        list_m[[i]] <- list_m[[i]][x, ]
    
    return(list_m)
}