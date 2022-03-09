# Reorder the color of a ggplot (in case of a missing modality)
order_color <- function(df, p, matched = NULL, collapse = FALSE,
                        colors = NULL) {
  J <- names(df)

  if (collapse) {
    n <- J
  } else {
    n <- J[-length(J)]
  }

  f <- ifelse(is.null(matched), scale_color_manual, scale_fill_manual)
  if (is.null(matched)) matched <- seq(n)

  return(p + f(
    values = color_group(J, colors)[matched],
    limits = n[matched],
    labels = n[matched]
  ))
}
