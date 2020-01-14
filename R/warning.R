warning <- function(message, call = sys.call(-1)) {
    base::warning(message, call. = FALSE, immediate. = TRUE)
}
