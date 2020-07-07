# p: a ggplot function
# hovertext : a boolean for the use of hovertext (if TRUE) as the attribute 
# to parse the onMouseOver text ("text" attribute, if FALSE)
# type: class of graphic among regular, var1D, perm
# p_perm: permutation plot object
modify_hovertext <- function(p, hovertext = TRUE, type = "regular", p_perm = NULL) {

    attr <- ifelse(hovertext, "hovertext", "text")
    # identify the order / id of the traces which corresponds to x- and y-axis 
    # (should be before the splitting function)

    if (type == "perm")
        k <- 2
    else if (type == "var1D")
        k <- 1
    else
        k <- seq(2)

    if (type == "perm") {
        n <- 2

    } else if (type == "regular") {

    traces <- which(lapply(
        p$x$data,
        function(x) length(grep("intercept", x$text)) == 1) == TRUE)

    # length of groups of points without traces and circle points
    n <- which(sapply(p$x$data, function(x)
        match("xintercept: 0", x$text) == 1)) - 1

    }else
        n <- length(p$x$data)

    for (i in seq(n)) {
        # For each lines of each group in the legend
        for (j in seq(length(p$x$data[[i]][attr][[1]]))) {

            # Distinguish each duplicate by splitting with "<br>" and separe
            # them in key/value by splitting with ": " (like a dictionnary)
            l_text <- lapply(
                as.list(strsplit(p$x$data[[i]][attr][[1]][j], "<br />")[[1]]),
                function(x) strsplit(x, ": ")[[1]])

            # keep only the (x, y) coordinates with the key df[, ]
            # and the response if exists
            l_text = unlist(lapply(l_text, function(x, y) {

                if ((type == "boot" && !is.character2(x[2]))
                    || x[1] %in% paste0("df[, ", k, "]"))
                        round(as.numeric(x[2]), 3)
                else if (type != "var1D" && x[1] == "resp")
                    x[2]
            }))

            if (type %in% c("regular", "boot")) {
                # do not print names because text = FALSE plots have'nt names
                name = ifelse(hovertext,
                            paste0("name: ", p$x$data[[i]]$text[j], "<br />"),
                            "")
                # Overwrite the onMouseOver text with the (x, y) coordinates and
                # the response if exists
                p$x$data[[i]][attr][[1]][j] <-
                    paste0(name,"x: ",l_text[1], "<br />y: ", l_text[2],
                                ifelse(
                                    length(l_text) == 3,
                                    paste0("<br />response: ", l_text[3]) ,
                                    ""
                                ))
            } else {
                if (type == "perm") {
                    res <- as.list(round(attributes(p_perm)$penalties[j, ], 3))
                    l_text <- paste0(
                        gsub("</?i>| ", "", p$x$layout$yaxis$title$text),
                        ": ",
                        l_text)
                    for (b in seq(length(names(res))))
                        l_text <- paste0(l_text, "<br />", paste0(names(res)[b], ": ", res[b]))
                }

                p$x$data[[i]][attr][[1]][j] <- l_text
            }

        }
    }

    # Remove the x- and y- axis onOverMouse
    if (type == "regular")
        (plotly::style(p, hoverinfo = "none", traces = traces))
    else
        p
}
