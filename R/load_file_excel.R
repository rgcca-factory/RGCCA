#' Creates a data frame from an Excel file loading
#
#' @param file A character giving the file name
#' @param sheet A character giving the sheet name
#' @param rownames An integer corresponding to the column number of the row
#' names (NULL otherwise)
#' @param header A bolean giving the presence or the absence of the header
#' @param num A bolean giving the presence or the absence of numerical values
#' @return A matrix containing the loaded file
# @examples
# \dontrun{
# load_file_excel("data/blocks.xlsx", "industry")
# }
# @export loadExcel
load_file_excel = function(
    file, 
    sheet = 1, 
    rownames = 1,
    header = TRUE,
    num = TRUE) {

    load_libraries("openxlsx")
    if (!is.null(rownames) && rownames < 1)
        rownames <- NULL

    df <- read.xlsx(
            file,
            sheet,
            colNames = header,
            na.strings = "NA")

    if (!is.null(rownames)) {
        names <- df[, rownames]
        df <- df[, -rownames, drop = FALSE]
    }
    
    if (num) 
        df <- as.data.frame(lapply(df, function(x) as.numeric(as.vector(x))))
    
    df <- as.matrix(df)
    
    if (!is.null(rownames)) row.names(df) <- names
    
    return(df)
}
