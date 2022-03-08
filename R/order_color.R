# Reorder the color of a ggplot (in case of a missing modality)
order_color <- function(df, p, matched = NULL, collapse = FALSE, colors = NULL) {
  J <- names(df)

  if (collapse) {
    n <- J
  } else {
    n <- J[-length(J)]
  }

  if (is.null(matched)) {
    matched <- seq(n)
    f <- "color"
  } else {
    f <- "fill"
  }

  func <- quote(get(paste("scale", f, "manual", sep = "_"))(
    values = color_group(J, colors)[matched],
    limits = n[matched],
    labels = n[matched]
  ))

  return(p + eval(as.call(func)))
}
