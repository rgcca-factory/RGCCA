#' Groups of color
#' 
#' Returns a color vector of equal size to the input vector
#' @param x A vector
#' @return A color vector of equal size to the input vector
#' @examples
#' color_group(seq(10))
#' @export
color_group <- function(x) {
    palette <-
        rep(
            c(
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
            ),
            10
        )
    palette[0:length(levels(as.factor(x)))]
}
