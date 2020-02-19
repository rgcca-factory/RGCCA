load_connection <- function(file = NULL, separator = "\t"){
    if (!is.null(file))
        load_file(
            file = file,
            separator = separator,
            rownames = NULL,
            header = FALSE
        )
}
