check_file <- function(f) {
    # Check the existence of a path f: A character giving the path of a file

    if (!file.exists(f))
        stop(paste(f, "does not exist."), exit_code = 101)

}
