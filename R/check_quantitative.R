#' Check if a dataframe contains no qualitative variables
#'
#' @param df A dataframe or a matrix
#' @param fo A character giving the name of the tested file
#' @param header A bolean giving either the presence (TRUE) or absence (FALSE) of
#'  a header
#' @param warn_separator A bolean to print warning for bad separator use
check_quantitative <- function(df, fo, header = FALSE, warn_separator = FALSE) {
    qualitative <- is.character2(df, warn_separator = TRUE)

    if (qualitative) {
        msg <- paste(
            fo,
            "contains qualitative data. Please, transform them in a disjunctive table."
        )

        if (!header)
            msg <- paste0(msg, "Possible mistake: header parameter is disabled, check if the file doesn't have one.")

        stop(paste(msg, "\n"), exit_code = 100)
    }

}

