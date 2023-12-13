#' The function subset_rows extracts rows from an object (vector, matrix, array,
#' data.frame).
#' @param x An object from which we want to extract rows
#' @param rows A set of rows
#' @noRd
subset_rows <- function(x, rows) {
  is.x.data.frame <- is.data.frame(x)
  if (is.x.data.frame) {
    row.names <- attr(x, "row.names")[rows]
  }
  if (is.vector(x)) {
    x <- x[rows]
  } else {
    x <- apply(x, -1, "[", rows)
  }
  if (is.x.data.frame) {
    x <- data.frame(x, row.names = row.names)
  }
  return(x)
}
