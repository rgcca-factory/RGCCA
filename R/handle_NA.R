#' Set up strategy to handle missing values
#'
#' @noRd
handle_NA <- function(blocks, NA_method = "nipals") {
  available_NA_methods <- c("complete", "nipals")
  na.rm <- FALSE
  if (!NA_method %in% available_NA_methods) {
    stop_rgcca(paste0(
      "NA_method ", NA_method, " is not implemented to handle missing values.",
      "Please select one among (",
      paste(available_NA_methods, collapse = ", "), ")."
    ))
  }
  if (NA_method == "complete") blocks <- intersection_list(blocks)
  if (NA_method == "nipals") {
    blocks <- blocks
    na.rm <- Reduce("||", lapply(blocks, function(x) any(is.na(x))))
  }
  return(list(blocks = blocks, na.rm = na.rm))
}
