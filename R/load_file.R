load_file <- function(
    file,
    file_text = file,
    separator = "\t",
    sheet = 1,
    rownames = 1,
    header = TRUE,
    one_column = FALSE,
    decimal = ".",
    rm_dup_rows = TRUE
    ) {

    match.arg(separator, c("\t", ",", ";"))
    match.arg(decimal, c(",", "."))
    check_boolean(header)
    if (!is.null(rownames))
        check_integer("rownames", rownames)
     if (is.character2(sheet))
        sheet <- paste0(sheet, collapse = " ")

    isXls <- length(grep("xlsx?", file))

    if (!isXls)
        file <- file_text

    check_file(file)
    check_size_file(file)

    # TODO: add automatic separator setting

    if (!isXls)
        load_file_text(file, separator, rownames, header, one_column, decimal, rm_dup_rows)
    else
        load_file_excel(file, sheet, rownames, header = header)

}
