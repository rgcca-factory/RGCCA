#' Creates a matrix from loading a file
#'
#' @inheritParams load_blocks
#' @param one_column A boolean for a file with  an only-one column
#' @return A matrix containing the loaded file
#' @examples
#' \dontrun{
#' load_file_text('data/agriculture.tsv')
#' }
load_file_text <- function(file, sep = "\t", rownames = 1, header = TRUE, one_column = FALSE) {

    if (!is.null(rownames) && rownames < 1)
        rownames <- NULL

    func <- function(x = rownames)
        as.matrix(read.table(
            file,
            sep = sep,
            header = header,
            row.names = x,
            na.strings = "NA",
            dec = ","
        ))

    tryCatch(
        f <- func(),
    error = function(e) {
        if (e$message == "duplicate 'row.names' are not allowed")
            f <<- func(NULL)
    })

    if (!one_column && NCOL(f) == 0)
        stop(paste(basename(file), "has an only-column. Check the separator."),
        exit_code = 102)

    return(f)
}
