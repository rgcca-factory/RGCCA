#' Set up strategy to handle missing values
#'
#' @noRd
handle_NA <- function(blocks, NA_method = "nipals") {
  na.rm <- FALSE
  if (NA_method == "complete") blocks <- intersection_list(blocks)
  if (NA_method == "nipals") {
    blocks <- blocks
    na.rm <- Reduce("||", lapply(blocks, function(x) any(is.na(x))))
  }
  return(list(blocks = blocks, na.rm = na.rm))
}
