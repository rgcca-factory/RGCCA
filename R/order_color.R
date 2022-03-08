# Reorder the color of a ggplot (in case of a missing modality)
order_color <- function(df, p, matched = NULL, collapse = FALSE,
                        colors = NULL) {
  J <- names(df)

  if (collapse) {
    n <- J
  } else {
    n <- J[-length(J)]
  }

  f <- ifelse(is.null(matched), "color", "fill")
  if (is.null(matched)) matched <- seq(n)

  func <- quote(get(paste("scale", f, "manual", sep = "_"))(
    values = color_group(J, colors)[matched],
    limits = n[matched],
    labels = n[matched]
  ))

  return(p + eval(as.call(func)))
}
