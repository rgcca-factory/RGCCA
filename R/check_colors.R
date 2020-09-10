check_colors <- function(colors, n = 3){

    if (!is.null(colors)){
        colors <- as.vector(colors)
        for (i in colors) {
            if (!is.na(i) && !(i  %in%  colors()) && is.character2(i)
                    && regexpr("^#{1}[a-zA-Z0-9]{6,8}$", i) < 1)
               stop_rgcca("colors must be in colors() or a rgb character.")
        }
    }
}
