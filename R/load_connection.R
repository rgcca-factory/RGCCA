load_connection <- function(file = NULL,
                            separator = "\t",
                            rownames = 1) {
  if (!is.null(file)) {
    load_file(
      file = file,
      separator = separator,
      rownames = rownames,
      header = TRUE
    )
  }
}
