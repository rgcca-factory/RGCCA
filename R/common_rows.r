# Keep only the rows with the same names among a list of dataframe
#
# @param list_m A list of matrix
# @return A list of matrix
common_rows <- function(list_m) {
  x <- Reduce(intersect, lapply(list_m, row.names))
  colnames_list <- lapply(list_m, colnames)
  for (i in seq(length(list_m)))
  {
    list_m[[i]] <- as.matrix(list_m[[i]][x, ])
    colnames(list_m[[i]]) <- colnames_list[[i]]
  }
  return(list_m)
}
