#' Plot permuation in 2D
#' 
#' Plot permuation in 2D
#' 
#' @inheritParams plot_var_2D
#' @param perm A permutation object from a RGCCA analyse
#' @param type An string giving the type of the index to look at (among 'crit' for
#'  the RGCCA criterion and 'zstat' for the pseudo Z-score)
#' @examples
#' data("Russett")
#' A = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
#'     politic = Russett[, 6:11] )
#' perm <- rgcca_permutation(A, nperm = 2, n_cores = 1)
#' plot_permut_2D(perm)
#' perm <- rgcca_permutation(A, p_c1 = TRUE, nperm = 2, n_cores = 1)
#' plot_permut_2D(perm)
#' @export
plot_permut_2D <- function(
    perm, 
    type = "zstat",
    cex = 1,
    cex_sub = 16 * cex,
    cex_point = 3 * cex,
    cex_lab = 19 * cex) {

    match.arg(type, c("crit", "zstat"))
    
    switch(
        type,
        "zstat" =  y_title <- "Z-score",
        "crit" = y_title <- "RGCCA criterion"
    )

    y <- unlist(perm[type])
    best <- which.max(y)
    n <- seq(nrow(perm$penalties))

    df <- as.data.frame(cbind(seq(NROW(perm$penalties)), y))
    rownames(df) <- sapply(n, function(x) paste(round(perm$penalties[x, ],2), collapse = ","))
    colnames(df) <- c("iter", type)

    axis <- function(margin){
        element_text(
            face = "italic",
            size = cex_lab * 0.75,
            margin = margin
        )
    }

    p <- ggplot(data = df, mapping = aes(x = df[, 1], y = df[, 2], ymin = 0)) + 
        theme_classic() +
        geom_line(size = 1) +
        labs(
            title = paste0(
                "Permutation scores \n(best value : ",
                paste(round(perm$penalties[best,], 2), collapse = ", "),
                ")"
            ), 
            x = "Index of combination",
            y = y_title
        ) +
        theme_perso(cex, cex_sub) +
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
                y = y[best],
                color = "red",
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
    else
        for (i in NCOL(perm$permcrit)) {
            p <- p + geom_line(
                aes(
                    x = df[, 1],
                    y = perm$permcrit[, i]
                ),
                col = "grey",
                size = 1
            )
        }

    return(p)

}
