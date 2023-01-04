#' Format bootstrap results
#'
#' Raw bootstrap results are stored in a list of n_boot elements where each
#' element of this list is a list of J matrices. Each of these matrices has
#' ncomp[j] columns where ncomp[j] is the number of components for block j.
#' The aim of this function is to reorganize this list into a data.frame.
#' @param W raw bootstrap results
#' @return A data.frame containing the results.
#' @noRd

format_bootstrap_list <- function(W, rgcca_res) {
  # Repeat values for variables, blocks and components
  grid <- do.call(rbind, lapply(seq_along(rgcca_res$a), function(j) {
    a <- rgcca_res$a[[j]]
    expand.grid(list(
      var = rownames(a),
      comp = seq_len(NCOL(a)),
      block = names(rgcca_res$a)[j]
    ), stringsAsFactors = FALSE)
  }))

  # Repeat values by adding the type and the number of the bootstrap sample
  grid <- merge(grid, expand.grid(list(
    type = c("weights", "loadings"),
    boot = seq_along(W)
  ), stringsAsFactors = FALSE), by = NULL)

  # Unlist the values into a new column of the grid
  grid$value <- unlist(W, use.names = FALSE)

  return(grid)
}
