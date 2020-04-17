
#' @importFrom rlang error_cnd

stop <- function(
    message = "",
    exit_code = 1) {

    base::stop(error_cnd(.subclass = exit_code, message = message))
}
        