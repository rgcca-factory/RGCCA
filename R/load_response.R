#' Create a matrix corresponding to the response
#'
#' @inheritParams load_blocks
#' @inheritParams set_connection
#' @return A matrix corresponding to the response
#' @examples
#' \dontrun{
#' blocks = lapply(seq(3), function(x) matrix(runif(47 * 5), 47, 5))
#' load_response (blocks, 'inst/extdata/response3.tsv')
#' }
#' @export
load_response <- function(
    blocks = NULL,
    file = NULL,
    sep = "\t",
    header = TRUE,
    rownames = 1) {

    response <- NULL

    if (!is.null(file)) {

        response <- load_file(
            file,
            sep = sep,
            rownames = rownames,
            header = header,
            one_column = TRUE
        )
    }

    check_response(response, blocks)

}
