subset_block_rows <- function(x, rows, drop = TRUE) {
  UseMethod("subset_block_rows")
}

#' @export
subset_block_rows.numeric <- function(x, rows, drop = TRUE) {
  return(x[rows])
}

#' @export
subset_block_rows.data.frame <- function(x, rows, drop = TRUE) {
  row.names <- attr(x, "row.names")[rows]
  x <- apply(x, -1, "[", rows, drop = drop)
  data.frame(x, row.names = row.names)
}

#' @export
subset_block_rows.array <- function(x, rows, drop = TRUE) {
  apply(x, -1, "[", rows, drop = drop)
}
