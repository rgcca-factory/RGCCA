#' Keep only the rows with the same names among a list of dataframe
#'
#' @param list_m A list of dataframe
#' @return A list of dataframe
keep_common_rows <- function(list_m) {

    names <- names(list_m)
    common_row <- row.names(list_m[[1]])

    for (i in 2:length(list_m))
        common_row <- common_row[common_row %in% row.names(list_m[[i]])]

    list_m <- lapply(
        seq(length(list_m)),
        function(x)
            list_m[[x]] <- list_m[[x]][common_row, ])

    names(list_m) <- names
    return(list_m)
}
