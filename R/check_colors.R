check_colors <- function(colors){
    
    if (is.null(colors)) {
        colors[1] <- "#A50026"
        colors[2] <- NA
        colors[3] <- "#313695"
    }else{
        if (length(colors) < 3)
            stop("color parameter must be at least of length 3.")
    }

    return(colors)
}