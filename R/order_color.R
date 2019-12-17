# Reorder the color of a ggplot (in case of a missing modality)
order_color <- function(df, p, matched = NULL, collapse = FALSE) {

    J <- names(df)
    if (is.null(matched)) {
        matched <- seq(length(J))
        f <- "color"
    } else
        f <- "fill"

    if (collapse)
        n <- J
    else
        n <- J[-length(J)]

    func <- quote(get(paste("scale", f, "manual", sep = "_"))(
        values = color_group(J)[matched],
        limits = n[matched],
        labels = n[matched]
    ))
    
    return(p + eval(as.call(func)))
}
