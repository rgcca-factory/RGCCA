# Retunr the block of each variable
#
# Get a vector of block names for each corresponding variable. The last block
# is considered as the superblock and ignored.
#
# @inheritParams plot2D
# @param df A list of matrix where their names are those of the blocks and the
# superblock and their rows are named after their variables
# @return A vector of character giving block names for each corresponding
# variable.
# @seealso \code{\link[RGCCA]{rgccad}}, \code{\link[RGCCA]{sgcca}}

get_bloc_var <- function(df, collapse = FALSE) {
  if (!collapse) {
    bl_names <- names(df)[-length(df)]
  } else {
    bl_names <- names(df)
  }

  res <- rep(
    bl_names,
    sapply(
      df[seq(length(df) - as.integer(!collapse))],
      function(x) NROW(as.matrix(x))
    )
  )

  names(res) <- unlist(lapply(df[bl_names], row.names))

  return(res)
}
