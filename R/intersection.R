#' Keeps subjects without missing values
#' @inheritParams rgccad
#' @return A list of blocks without missing values and having the same row names
#' @title intersection_list
#' @noRd
intersection_list <- function(A) {
  # Find rows without missing values in each block
  valid_rows <- lapply(A, complete.cases)
  # Take the intersection
  common_valid_rows <- apply(
    matrix(unlist(valid_rows), length(valid_rows[[1]]), length(valid_rows)),
    1, prod
  )
  if (sum(common_valid_rows) <= 3) {
    stop_rgcca(paste0(
      "Less than 3 subjects have no missing values, choose",
      " another missing value handling method or work on ",
      "your dataset."
    ))
  }
  # Extract the rows from the different blocks
  lapply(A, function(x) subset_rows(x, as.logical(common_valid_rows)))
}
