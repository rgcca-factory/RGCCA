#' Histogram of Average Variance Explained
#'
#' Histogram of the model quality (based on Average Variance Explained)
#' for each blocks and sorted in decreasing order
#'
#' @inheritParams plot_ind
#' @inheritParams plot_histogram
#' @seealso \code{\link[RGCCA]{rgcca}}, \code{\link[RGCCA]{sgcca}}
#' @examples
#' random_val = function(y=1) lapply(seq(4),
#' function(x) matrix(runif(4), y, 2))
#' rgcca_out = list(AVE = list(AVE_X = random_val()),
#'      a = random_val(2), ncomp = rep(2, 4))
#' names(rgcca_out$a) <- LETTERS[seq(4)]
#' library("ggplot2")
#' for(i in seq(1,4))
#' names(rgcca_out$AVE$AVE_X[[i]]) <- c(1,2)
#' plot_ave(rgcca_out)
#' @export
plot_ave <- function(
    rgcca,
    cex = 1,
    cex_sub = 16 * cex,
    cex_axis = 10 * cex) {

    if (is(rgcca, "pca")) {
        rgcca$AVE$AVE_X = rgcca$AVE$AVE_X[1]
        rgcca$ncomp = rgcca$ncomp[1]
        rgcca$a = rgcca$a[1]
    }
    
    names(rgcca$AVE$AVE_X) <- NULL
    ave <- 100 * unlist(rgcca$AVE$AVE_X)
    blocks <- factor(unlist(lapply(seq(length(names(rgcca$a))),
            function(x) rep(names(rgcca$a)[x], rgcca$ncomp[x]))),
        levels = names(rgcca$a))
    ncomp <- as.factor(names(ave))

    y_ave_cum <- lapply(
        lapply(rgcca$AVE$AVE_X, 
            function(x) round(100 * cumsum(x), 1)), 
        function(x) c(0, x))
    y_ave_cum <- unlist(lapply(y_ave_cum, function(x)
            unlist(lapply(seq(length(x)),
                function(i) (x[i - 1] + x[i]) / 2))))

    ave_label <- unlist(lapply(rgcca$AVE$AVE_X, function(x)
            round(100 * x, 1)))
    ave_label[ave_label < max(y_ave_cum) / 20] <- ""

    df <- data.frame(ave, blocks, ncomp, stringsAsFactors = FALSE)

    p <- ggplot(data = df, 
        aes(
            x = blocks,
            y = ave,
            fill = ncomp,
            label = ave_label
        ))

    p <- plot_histogram(
        p, 
        df, 
        "Average Variance Explained",
        cex = cex,
        cex_sub = cex_sub,
        cex_axis = cex_axis) +
    scale_fill_manual(
        values = color_group(levels(df$ncomp)),
        labels = gsub("comp", " ", levels(df$ncomp))) + 
    geom_col(position = position_stack(reverse = TRUE)) +
    labs(subtitle = print_comp(rgcca, outer = TRUE)) +
    geom_text(aes(y = y_ave_cum),  cex = 3.5 * cex, color = "white") +
    labs(fill = "Components")

    return(p)
}
