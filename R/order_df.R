# Ordered rows of a data.frame in decreasing order.
#
# @inheritParams plot_var_1D
# @param df A dataframe
# @param allCol A logical indicating to keep all the columns in the data.frame
# @return A data.frame with ordered values

order_df <- function(df, comp = 1, allCol = TRUE) {
  ordered <- order(abs(df[, comp]), decreasing = TRUE)

  if (allCol) comp <- seq(NCOL(df))
  res <- df[ordered, comp]

  if (!allCol) {
    names(res) <- row.names(df)[ordered]
  }

  return(res)
}
