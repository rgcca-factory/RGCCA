 #' @importFrom rlang error_cnd
 stop_rgcca <- function(
     message = "",
     exit_code = 1) {
    # load_libraries("rlang")
     base::stop(error_cnd(.subclass = exit_code, message = message,call.=T,domain=NA))
    # stop(message,call.=FALSE,domain=NA)
 }
#'         