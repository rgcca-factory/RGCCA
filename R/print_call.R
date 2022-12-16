print_call <- function(call) {
  ### Print parameters of the function
  cat("Call: ")
  names_call <- c(
    "method", "superblock", "scale", "scale_block", "init",
    "bias", "tol", "NA_method", "ncomp", "response"
  )

  char_to_print <- vapply(names_call, function(name) {
    if (name == "ncomp") {
      value <- (paste(call$ncomp, sep = "", collapse = ","))
      value <- paste0("c(", value, ")")
    } else {
      value <- call[[name]]
    }
    quo <- ifelse(is.character(value) && (name != "ncomp"), "'", "")
    if (is.null(value)) value <- "NULL"
    paste0(name, "=", quo, value, quo)
  }, FUN.VALUE = character(1))
  cat(paste(char_to_print, collapse = ", "), "\n")

  ### Print number of blocks
  cat("There are J =", NCOL(call$connection), "blocks.", fill = TRUE)

  ### Print design matrix
  cat("The design matrix is:\n")
  print(call$connection)

  ### Print scheme
  cat("\n")
  if (is.function(call$scheme)) {
    cat("The", deparse(call$scheme), "scheme is used.", fill = TRUE)
  } else {
    cat("The", call$scheme, "scheme is used.", fill = TRUE)
  }
}
