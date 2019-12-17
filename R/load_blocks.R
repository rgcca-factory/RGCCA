#' Create a list of matrix from loading files corresponding to blocks
#'
#' @param file A character giving the path of a file
#' @param names A character giving a list of names for the blocks
#' @param sep A character giving the column separator
#' @param header A bolean giving the presence or the absence of the header
#' @param rownames An integer corresponding to the column number of the row
#' names (NULL otherwise)
#' @return A list matrix corresponding to the blocks
#' @examples
#' \dontrun{
#' load_blocks (TRUE,
#'     "inst/extdata/agriculture.tsv,inst/extdata/industry.tsv,inst/extdata/politic.tsv",
#'     "agric,ind,polit")
#' }
#' @export
load_blocks <- function(file,
    names = NULL,
    sep = "\t",
    header = TRUE,
    rownames = 1) {

    # Parse args containing files path
    isXls <- (length(grep("xlsx?", file)) == 1)
    # test if extension filename is xls
    if (!isXls)
    # if it is not, parse the name of file from the arg list
        block_filenames <- char_to_list(file)
    else {
        # # if xls, check file exists
        # check_file(file)
        # # load the xls
        # wb = loadWorkbook(file)
        # # load the blocks
        # block_filenames = names(getSheets(wb))
    }

    # Parse optional names of blocks
    if (!is.null(names)) {
        # default name is filename, otherwise, the user could name the blocs
        block_names <- char_to_list(names)
        invisible(
            check_size_blocks(
                block_filenames,
                "names",
                block_names
        ))
    }

    # Load each dataset
    blocks <- list()
    for (i in seq(length(block_filenames))) {
        if (!isXls) {
            fi <- block_filenames[i]
            check_file(fi)
        }

        #Get names of blocs
        if (!is.null(names))
            # names of blocks are those parsed from args
            fo <- get_filename(block_names[i])
        else {
            if (!isXls)
                # if not xls, the name is the files without the extension .tsv
                fo <- get_filename(fi)
            else
                # for xls, the names are those of the sheets
                fo <- block_filenames[i]
        }

        df <- load_file(file, fi, sep, block_filenames[i], rownames, header)

        check_quantitative(df, fo, header)
        blocks[[fo]] <- df
    }

    blocks <- check_blocks(blocks, init = TRUE)

}
