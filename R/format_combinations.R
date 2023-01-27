#' Small utility function to format a matrix of parameters
#' @noRd
format_combinations <- function(par_value) {
  combinations <- apply(
    format(par_value, digits = 2), 1, paste0, collapse = "/"
  )
  # If parameters are too long, there are replaced with "Set x"
  # The same is done if rounding to 2 digits leads to the same values
  to_set <- (nchar(combinations[1]) > 15) |
    (length(unique(combinations)) < NROW(par_value))
  if (to_set) {
    combinations <- paste("Set ", sep = "", seq_len(NROW(par_value)))
  }
  return(combinations)
}
