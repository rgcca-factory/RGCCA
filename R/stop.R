stop <- function(
    message = "",
    exit_code = 1) {

    base::stop(rlang::error_cnd(.subclass = exit_code, message = message))
}
        