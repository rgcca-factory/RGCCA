#' Plot permuation in 2D
#' 
#' Plot permuation in 2D
#' 
#' @inheritParams plot_var_2D
#' @param perm A permutation object from a RGCCA analyse
#' @param bars Among "points", "stderr" or "sd": representation of the variability
#' @param type An string giving the type of the index to look at (among 'crit' for
#'  the RGCCA criterion and 'zstat' for the pseudo Z-score)
#' @param cex general size of the text
#' @param cex_main = 25 * cex, size of the main text (title)
#' @param cex_sub = 16 * cex, size of the subtitle text 
#' @param cex_point = 3 * cex, size of the point
#' @param cex_lab = 19 * cex, size of the labels
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 geom_vline

plot_permut_2D <- function(
    perm, 
    type = "crit",
    cex = 1,
    title = NULL,
    cex_main = 14 * cex,
    cex_sub = 12 * cex,
    cex_point = 3 * cex,
    cex_lab = 10 * cex,
    bars="points"
    ) {

    match.arg(type, c("crit", "zstat"))
    for (i in c("cex", "cex_main", "cex_sub", "cex_point", "cex_lab"))
        check_integer(i, get(i))
    
    switch(
        type,
        "zstat" =  y_title <- "Z-score",
        "crit" = y_title <- "RGCCA criterion"
    )

    check_ncol(list(perm$zstat), 1)

    y <- unlist(perm[type])
    best <- which.max(unlist(perm["zstat"]))
    y_best <- y[best]
    n <- seq(nrow(perm$penalties))

    df <- as.data.frame(cbind(seq(NROW(perm$penalties)), y))
    rownames(df) <- n
    colnames(df) <- c("iter", type)

    axis <- function(margin){
        element_text(
            face = "italic",
            size = cex_lab * 0.75,
            margin = margin
        )
    }

    if (is.null(title))
        title <- paste0(
                "Permutation scores \n(best value : ",
                paste(round(perm$penalties[best,], 2), collapse = ", "),
                ")"
            )
    else
        title <- paste0(title, collapse = " ")

    p <- ggplot(data = df, mapping = aes(x = df[, 1], y = df[, 2], ymin = 0)) + 
        theme_classic() +
        geom_line(size = 1) +
        labs(
            title = title, 
            x = "Combinations",
            y = y_title
        ) +
        theme_perso(cex, cex_main, cex_sub) +
        theme(
            axis.text = element_text(size = 10, face = "bold"),
            axis.title.y = axis(margin(0, 20, 0, 0)),
            axis.title.x = axis(margin(20, 0, 0, 0)),
            axis.line = element_line(size = 1),
            axis.ticks  = element_line(size = 1),
            axis.ticks.length = unit(2, "mm"),
            legend.position = "none"
        ) + 
        geom_point(
            mapping = aes(
                x = best,
                y = y_best,
                color = I("red"),
                shape = I(3)
            ),
            size = 5
        ) +
        geom_vline(
            size = 1,
            color = "red",
            xintercept = best
        )

    if (type == "zstat")
        p <- p + geom_hline(
            size = 1,
            color = "grey",
            linetype = "dashed",
            yintercept = c(1.96, 2.58, 3.29)
        )
    else {
        dft <- NULL
        for (i in seq(NCOL(perm$permcrit))) {
            x <- df[, 1]
            y <- perm$permcrit[, i]
            dfi <- cbind(x,y)
            dft <- rbind(dft,dfi)
        }
        dft <- as.data.frame(dft)
     
        if (bars == "points")
            p <- p + geom_point(data = dft,aes(x = dft[,1], y = dft[,2]), colour = "green", size = 0.8)
        if (bars == "sd") {
        }
        if (bars == "stderr") {
        }
        if (bars == "ci") {
        }
        if (bars == "quantile") {
        }

    }

    p <- p + scale_x_continuous(breaks = 1:nrow(df), labels = rownames(df))
    attributes(p)$penalties <- perm$penalties

    return(p)
}
