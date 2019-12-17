# Print warning if file size over
check_size_file <- function(filename) {
    size <- file.size(filename)
    if (size > 5e+06)
        # warning(paste0('The size of ', filename, ' is over 5 Mo (',
        #  round(size / 1E6, 1), ' Mo). File loading could take some times...'),
        message("File loading in progress ...")
}
