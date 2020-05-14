 #' @importFrom rlang error_cnd
 stop_rgcca <- function(
     message = "",
     exit_code = 1) {
    # load_libraries("rlang")
     base::stop(error_cnd(.subclass = exit_code, message = paste0("\n",message),call.=T,domain=NA))
    # stop(message,call.=FALSE,domain=NA)
 }
#'         