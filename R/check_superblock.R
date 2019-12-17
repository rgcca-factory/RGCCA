check_superblock <- function(is_supervised = NULL, is_superblock = NULL) {
    if (!is.null(is_supervised)) {
        warn_connection("supervized method with a response")
        if (is_superblock) {
            if (!is.null(is_superblock))
                warning("In a supervised mode, the superblock corresponds to the response.")
        }
        return(FALSE)
    }else
        return(isTRUE(is_superblock))
}
