#' Create a matrix corresponding to the response
#'
#' @inheritParams load_blocks
#' @inheritParams set_connection
#' @return A matrix corresponding to the response
#' @examples
#' \donttest{
#' blocks = lapply(seq(3), function(x) matrix(runif(47 * 5), 47, 5))
#' load_response (blocks, 'inst/extdata/response3.tsv')
#' }
#' @export
load_response <- function(
    blocks = NULL,
    file = NULL,
    separator = "\t",
    header = TRUE,
    rownames = 1,
    decimal = ".") {

    response <- NULL

    if (!is.null(file)) {

        response <- load_file(
            file,
            separator = separator,
            rownames = rownames,
            header = header,
            one_column = TRUE,
            decimal = decimal
        )
    }else
        check_blocks(blocks, , allow_unnames = FALSE, n = 1)

    check_response(response, blocks)

}
