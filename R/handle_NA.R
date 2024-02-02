#' Set up strategy to handle missing values
#'
#' @noRd
handle_NA <- function(blocks, NA_method = "na.ignore") {
  na.rm <- FALSE
  if (NA_method == "na.omit") blocks <- intersection_list(blocks)
  if (NA_method == "na.ignore") {
    blocks <- blocks
    na.rm <- any(is.na(unlist(blocks, use.names = FALSE)))
  }
  return(list(blocks = blocks, na.rm = na.rm))
}
