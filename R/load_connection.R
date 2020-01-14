load_connection <- function(file = NULL, separator = "\t"){
    if (!is.null(file))
        load_file(
            file = file,
            sep = separator,
            rownames = NULL,
            header = FALSE
        )
}
