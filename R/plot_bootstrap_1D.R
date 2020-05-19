#' Plot a bootstrap in 1D
#'
#' Histogram of the best variables of an RGCCA bootstrap with, on the x-axis,
#' the number of non-zero occurrences (SGCCA) or the bootstrap-ratio
#' (mean/sd; RCCA). The bars are colored according to the average weight of
#' the boostrap  (according to an ascending gradient from red to blue)
#' @inheritParams plot_histogram
#' @inheritParams get_bootstrap
#' @param df_b Result of get_bootstrap functions or dataframe #TODO
#' @param b A matrix of boostrap
#' @param x A character for the column used in the plot
#' @param y A character for the column to color the bars
#' @param n_mark An integer giving the number maximum of top variables
#' @param ... Other parameters (see plot_histogram)
#' @examples
#' library(RGCCA)
#' data("Russett")
#' blocks = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
#'     politic = Russett[, 6:11] )
#' rgcca_out = rgcca(blocks, sparsity = 0.75, type = "sgcca")
#' boot = bootstrap(rgcca_out, 2, n_cores = 1)
#' plot_bootstrap_1D(boot, n_cores = 1)
#' rgcca_out = rgcca(blocks)
#' boot = bootstrap(rgcca_out, 2, n_cores = 1)
#' selected.var = get_bootstrap(boot, n_cores = 1,display_order=TRUE)
#' plot_bootstrap_1D(boot, n_cores = 1)
#' plot_bootstrap_1D(df_b = selected.var)
#' @export
#' @importFrom ggplot2 ggplot
plot_bootstrap_1D <- function(
    b = NULL,
    df_b = NULL,
    x = "estimate",
    y = "occurrences",
    n_mark = 50,
    title = NULL, 
    colors = NULL,
    comp = 1,
    bars = "sd",
    i_block = length(b$bootstrap[[1]]),
    collapse = FALSE,
    n_cores = parallel::detectCores() - 1,
    ...) {

    if (missing(b) && missing(df_b))
        stop_rgcca("Please select a bootstrap object.")
    if (!is.null(b)) {
        df_b <- get_bootstrap(b, comp, i_block, collapse, n_cores, bars=bars,display_order = TRUE)
    }
    if (!is.null(df_b))
        stopifnot(is(df_b, "df_bootstrap"))
    check_integer("n_mark", n_mark)
    if (is.null(title))
        title <- paste0(attributes(df_b)$indexes[[x]],
                   "\n(",
                   attributes(df_b)$n_boot,
                   " bootstraps)")
    if (is.null(colors)) {
        if (!(y %in% c("occurrences", "sign")))
            colors <- c(color_group(seq(3))[1],  "gray", color_group(seq(3))[3])
        else
            colors <- c(color_group(seq(3))[1], color_group(seq(3))[3])
    }

    lower_band <- NULL -> upper_band
    check_ncol(list(df_b), 1)

    set_occ <- function(x) {
        match.arg(x, names(attributes(df_b)$indexes))
        if (x == "occurrences" && !x %in% colnames(df_b))
            return("sign")
        else
            return(x)
    }

    x <- set_occ(x)
    y <- set_occ(y)

    if (y == "sign") 
        group = seq(2)
    else
        group = NA

    df_b_head <- head(
        data.frame(
            order_df(df_b[, -NCOL(df_b)], x, allCol = TRUE),
            order = NROW(df_b):1),  n_mark)
    df_b_head<-df_b_head[df_b_head[,"sd"]!=0,]
    class(df_b_head) <- c(class(df_b), "d_boot1D")
    
    p <- ggplot(
        df_b_head,
        aes(x = order,
            y = df_b_head[, x],
            fill = df_b_head[, y]))

    p <- plot_histogram(
        p,
        df_b_head,
        title,
        group,
        colors,
        ...) +
    labs(fill = attributes(df_b)$indexes[[y]])

    if (x == "estimate")
        p <- p +
            geom_errorbar(aes(ymin = lower_band, ymax = upper_band,width=0.5))
    if(x =="occurrences")
    {
        n_boot=ifelse(!is.null(dim(b[[1]][[1]][[1]])),dim(b[[1]][[1]][[1]])[2],length(b[[1]][[1]][[1]]))
        nvar=length(b$bootstrap[[1]][[i_block]][,1])
        avg_n_occ=sum(df_b$occurrences)/n_boot
        probComp= avg_n_occ/nvar
        
        q1=qbinom(size=n_boot,prob=probComp,p=1-0.05/nvar)
        q2=qbinom(size=n_boot,prob=probComp,p=1-0.01/nvar)
        q3=qbinom(size=n_boot,prob=probComp,p=1-0.001/nvar)
        
        p <-p+geom_hline(yintercept = c(q1,q2,q3),col=c("red","black","green"))
        p
        
    }
    return(p)
}
