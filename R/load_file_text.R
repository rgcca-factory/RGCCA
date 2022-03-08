# Creates a matrix from loading a file
#
# @inheritParams load_blocks
# @param one_column A boolean for a file with  an only-one column
# @return A matrix containing the loaded file
# @examples
# \dontrun{
# load_file_text('data/agriculture.tsv')
# }
load_file_text <- function(file,
                           separator = "\t",
                           rownames = 1,
                           header = TRUE,
                           one_column = FALSE,
                           decimal = ".") {
  if (!is.null(rownames) && rownames < 1) {
    rownames <- NULL
  }

  func <- function(x = rownames) {
    as.matrix(read.table(
      file,
      sep = separator,
      header = header,
      row.names = x,
      na.strings = "NA",
      dec = decimal
    ))
  }

  tryCatch(
    f <- func(),
    error = function(e) {
      msg <- "duplicate 'row.names' are not allowed"
      if (e$message == msg) {
        message(paste0(msg, "; rownames have been removed from dataset."))
        f <<- func(NULL)
      }
    }
  )

  if (!one_column && NCOL(f) == 0) {
    stop_rgcca(paste(basename(file), "has an only-column. Check the separator."),
      exit_code = 102
    )
  }

  return(f)
}
