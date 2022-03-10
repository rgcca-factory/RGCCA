stop_rgcca <- function(...,
                       exit_code = "1",
                       call = NULL) {
  base::stop(
    structure(
      class = c(exit_code, "simpleError", "error", "condition"),
      list(message = paste0(...), call. = NULL)
    )
  )
}
