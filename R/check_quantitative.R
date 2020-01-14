#' Check if a dataframe contains no qualitative variables
#'
#' @param df A dataframe or a matrix
#' @param fo A character giving the name of the tested file
#' @param h A bolean giving either the presence (TRUE) or absence (FALSE) of
#'  a header
#' @examples
#' df = matrix(runif(20), 10, 2)
#' check_quantitative(df, 'data')
#' \dontrun{
#' df[,2] = LETTERS[seq(10)]
#' check_quantitative(df, 'data', TRUE)
#' # Error
#' }
#' @export
check_quantitative <- function(df, fo, h = FALSE) {
    qualitative <- is.character2(df)

    if (qualitative) {
        msg <- paste(
            fo,
            "contains qualitative data. Please, transform them in a disjunctive table."
        )

        if (!h)
            msg <- paste0(msg, "Possible mistake: header parameter is disabled, check if the file doesn't have one.")

        stop(paste(msg, "\n"), exit_code = 100)
    }

}

