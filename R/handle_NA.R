available_NA_methods <- c("complete", "nipals")

handle_NA <- function(blocks, method = "nipals") {
  if (!method %in% available_NA_methods) stop_rgcca(paste0(
    "method ", method, " is not implemented to handle missing values.",
    "Please select one among (", paste(available_NA_methods, collapse = ", "), ")."
  ))
  if (method == "complete") blocks = intersection_list(blocks)
  if (method == "nipals")   blocks = blocks
  return(blocks)
}
