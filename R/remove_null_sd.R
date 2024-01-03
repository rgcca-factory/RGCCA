#' Remove columns having a 0 standard deviation
#'
#' @param list_m A list of dataframe
#' @param column_sd_null Either NULL or a list of named vectors. If NULL, the
#' function will search for variables with null variance in each block. If not
#' NULL, this list defines for each block the index of the variables that are
#' of null variance (see the 'Value' section for more details about the content
#' of this list). In both cases, these variables are removed.
#' @return \item{list_m}{A list of dataframe.}
#' @return \item{column_sd_null}{Either NULL, if not a single variable was
#' removed, or a list of the same size as the number of blocks. In the last
#' situation, each element of this list is again NULL if not a single variable
#' was removed from the current block, or a named vector indicating the former
#' index of the removed variables along with their name.}
#' @noRd

remove_null_sd <- function(list_m, column_sd_null = NULL, groups = NULL) {
  names <- names(list_m)
  
  if (is.null(column_sd_null)) {
    column_sd_null <- lapply(
      list_m,
      function(x) {
        which(apply(x, 2, function(y) {
          if (mode(y) != "character") {
            if (!is.null(groups)){
              # in multi-group mode, variables with NA in a 
              # whole group but not in all groups are kept
              res <- all(sd(y[!is.na(y)]) == 0)
            } else {
              res <- all(is.na(y)) || (sd(y[!is.na(y)]) == 0)
            }
          } else {
            res <- FALSE
          }
          return(res)
        }))
      }
    )
  }
  
  # in multi-group mode, variables with null sd in at 
  # least one group are removed in all groups
  if (!is.null(groups)) {
    all_cols_sd_null <- unique(unlist(column_sd_null, use.names = FALSE))
    column_sd_null[seq_along(column_sd_null)] <- rep(list(all_cols_sd_null), 
                                                     times = length(column_sd_null))
  }
  
  blocks_index <- seq(1, length(list_m))[
    unlist(
      lapply(
        column_sd_null,
        function(x) length(x) > 0
      )
    )
  ]
  list_m <- lapply(
    seq_along(list_m),
    function(x) {
      if (x %in% blocks_index) {
        list_m[[x]][, -column_sd_null[[x]], drop = FALSE]
      } else {
        list_m[[x]]
      }
    }
  )
  
  names(list_m) <- names
  return(list(list_m = list_m, column_sd_null = column_sd_null))
}
