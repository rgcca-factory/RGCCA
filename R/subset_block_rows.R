subset_block_rows <- function(x, rows, drop = TRUE) {
  UseMethod("subset_block_rows")
}

#' @export
subset_block_rows.numeric <- function(x, rows, drop = TRUE) {
  return(x[rows, drop = drop])
}

#' @export
subset_block_rows.data.frame <- function(x, rows, drop = TRUE) {
  row.names <- attr(x, "row.names")[rows]
  x <- apply(x, -1, "[", rows, drop = drop)
  data.frame(x, row.names = row.names)
}

#' @export
subset_block_rows.array <- function(x, rows, drop = TRUE) {
  dim_x <- dim(x)
  dn <- dimnames(x)
  x <- apply(x, -1, "[", rows, drop = drop)
  if (!drop && length(dim(x)) < length(dim_x)) {
    x <- array(x, dim = c(1, dim_x[-1]))
    rownames(x) <- ifelse(is.numeric(rows), dn[[1]][rows], rows)
    dimnames(x)[-1] <- dn[-1]
  }
  return(x)
}
