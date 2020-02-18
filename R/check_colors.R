check_colors <- function(colors, n = 3){

    if (is.null(colors)) {
        colors <- "#A50026"
        colors[2] <- NA
        colors[3] <- "#313695"
    }else{
        colors <- as.vector(colors)
        for (i in colors) {
            if (!is.na(i) && !(i  %in%  colors()) && is.character2(i)
                    && regexpr("^#{1}[a-zA-Z0-9]{6,8}$", i) < 1)
               stop("colors must be in colors() or a rgb character.") 
        }
        # if (length(colors) < n)
        #     stop(paste0("color parameter must be at least of length ", n, "."))
    }

    return(colors)
}
