stop <- function(
    message,
    exit_code = "1",
    call = NULL) {

    base::stop(
        structure(
            class = c(exit_code, "simpleError", "error", "condition"),
            list(message = message, call. = NULL)
    ))
}
