#' Plot permutation in 2D
#' 
#' Plot permutation in 2D
#' 
#' @inheritParams plot_var_2D
#' @inheritParams plot2D
#' @inheritParams get_bootstrap
#' @param perm A permutation object (see \code{\link[RGCCA]{rgcca_permutation}})
#' @param bars A character giving the variability among "points" (all the points are shown), "sd" (standard deviations bars) , "stderr"(bars of standard deviation divided by sqrt(n)) or "quantile" (for the 0.05-0.95 quantiles bars)
#' @param type A character giving the type of the index to look at (among 'crit' for
#'  the RGCCA criterion and 'zstat' for the pseudo Z-score)
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 geom_vline 
#' @importFrom stats aggregate
plot_permut_2D <- function(
    perm, 
    type = "crit",
    cex = 1,
    title = NULL,
    cex_main = 14 * cex,
    cex_sub = 12 * cex,
    cex_point = 3 * cex,
    cex_lab = 10 * cex,
    bars = "points",
    colors = c("red", "grey")
    ) {

    xend <- yend <- NULL
    stopifnot(is(perm, "permutation"))
    match.arg(type, c("crit", "zstat"))
    match.arg(bars,c("points","sd","stderr","quantile"))
    for (i in c("cex", "cex_main", "cex_sub", "cex_point", "cex_lab"))
        check_integer(i, get(i))
    check_colors(colors)
    if (length(colors) < 2)
        colors <- rep(colors, 2)

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
                "Permutation scores (",perm$call$n_run, " runs) \n Best parameters : ",
                paste(round(perm$penalties[best,], 2), collapse = ", ")
            )
    else
        title <- paste0(title, collapse = " ")

    p <- ggplot(data = df, mapping = aes(x = df[, 1], y = df[, 2], ymin = 0)) + 
        theme_classic() +
        geom_line(size = 0.5) +
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
            axis.line = element_line(size = 0.5),
            axis.ticks  = element_line(size = 0.5),
            axis.ticks.length = unit(2, "mm"),
            legend.position = "none"
        )

    if (type == "zstat")
        p <- p + geom_hline(
            size = 0.5,
            color = colors[2],
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
            p <- p + geom_point(data = dft,aes(x = dft[,1], y = dft[,2]), colour = colors[2], size = 0.8)
         if (bars == "sd") {
             tab=aggregate(dft,by=list(dft[,1]),sd)
             tab2=aggregate(dft,by=list(dft[,1]),mean)
             dat=data.frame(x=tab[,1],y=tab2[,"y"]-tab[,"y"],xend=tab[,1],yend=tab2[,"y"]+tab[,"y"])
                p <- p+ geom_point(data=tab2,aes(x=tab2[,1],y=tab2[,3]),colour=colors[2])
               p <- p + geom_segment(data=dat,aes(x=x,y=y,xend=xend,yend=yend),colour=colors[2],size=0.5)
         }
         if(bars == "stderr")
         {
             tab=aggregate(dft,by=list(dft[,1]),function(x){return(sd(x)/sqrt(length(x)))})
             tab2=aggregate(dft,by=list(dft[,1]),mean)
             dat=data.frame(x=tab[,1],y=tab2[,"y"]-tab[,"y"],xend=tab[,1],yend=tab2[,"y"]+tab[,"y"])
             p <- p+ geom_point(data=tab2,aes(x=tab2[,1],y=tab2[,3]),colour=colors[2])
             p <- p + geom_segment(data=dat,aes(x=x,y=y,xend=xend,yend=yend),colour="green",size=0.5)
         }
        if(bars=="quantile")
        {
            tabq1=aggregate(dft,by=list(dft[,1]),function(x){return(quantile(x,0.05))})
            tabq2=aggregate(dft,by=list(dft[,1]),function(x){return(quantile(x,0.95))})
            tab2=aggregate(dft,by=list(dft[,1]),mean)
            dat=data.frame(x=tab2[,1],y=tabq1[,"y"],xend=tab2[,1],yend=tabq2[,"y"])
            p <- p+ geom_point(data=tab2,aes(x=tab2[,1],y=tab2[,3]),colour=colors[2])
            p <- p + geom_segment(data=dat,aes(x=x,y=y,xend=xend,yend=yend),colour=colors[2],size=0.5)
        }
   

    }

    p <- p + geom_line(data = df, mapping = aes(x = df[, 1], y = df[, 2])) +
        scale_x_continuous(breaks = 1:nrow(df), labels = rownames(df)) + 
        theme(plot.title = element_text(vjust=5), plot.margin = margin(5, 0, 0, 0, "mm")) + 
        geom_point(
            mapping = aes(
                x = best,
                y = y_best,
                color = I(colors[1]),
                shape = I(3)
            ),
            size = 5
        ) +
        geom_vline(
            size =  0.5,
            color = colors[1],
            xintercept = best
        )
    attributes(p)$penalties <- perm$penalties

    return(p)
}
