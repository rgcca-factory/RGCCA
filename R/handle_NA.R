#' Set up strategy to handle missing values
#'
#' @noRd
handle_NA <- function(blocks, NA_method = "na.ignore") {
  na.rm <- FALSE
  if (NA_method == "na.omit") blocks <- intersection_list(blocks)
  if (NA_method == "na.ignore") {
    blocks <- blocks
    na.rm <- Reduce("||", lapply(blocks, function(x) any(is.na(x))))
  }
  return(list(blocks = blocks, na.rm = na.rm))
}
