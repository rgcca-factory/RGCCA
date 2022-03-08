# Groups of color
#
# Returns a color vector of equal size to the input vector
# @inheritParams plot2D
# @param x A vector
# @return A color vector of equal size to the input vector

color_group <- function(x, colors = NULL) {
  if (is.null(colors)) {
    colors <- c(
      "#cd5b45",
      "#71ad65",
      "#3c78b4",
      "#ffc600",
      "#b448af",
      "#9d9d9d",
      "#abcaef",
      "#4a6f43",
      "#f0e500",
      "#efb8f0",
      "black",
      "#d6d6d6"
    )
  }

  colors <- rep(colors, 10)
  cols <- colors[seq(length(levels(as.factor(x))))]

  # names(cols)=levels(as.factor(x))
  return(cols)
}
