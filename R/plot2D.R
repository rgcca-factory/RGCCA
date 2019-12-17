#' Plot of components space
#'
#' Plots RGCCA components in a bi-dimensional space
#'
#' @inheritParams plot_ind
#' @inheritParams plot_var_2D
#' @param df A dataframe
#' @param title A character with the name of the space (either "Variables" or
#' "Samples")
#' @param group A vector of character with levels used to color the points
#' @param name_group A character giving the type of groups (either "Blocs" or
#' "Response")
#' @param p A ggplot object
#' @param colours A vectof of character to color quantitative dat
#' @examples
#' df = as.data.frame(matrix(runif(20*2, min = -1), 20, 2))
#' AVE = lapply(seq(4), function(x) runif(2))
#' rgcca_out = list(AVE = list(AVE_X = AVE))
#' plot2D(rgcca_out, df, "Samples", rep(c("a","b"), each=10), "Response")
#' @export
plot2D <- function(
    rgcca,
    df,
    title = "Biplot",
    group,
    name_group = "Response",
    compx = 1,
    compy = 2,
    i_block = NULL,
    p = NULL,
    text = TRUE,
    i_block_y = i_block,
    colours = c("blue", "gray", "#cd5b45"),
    collapse = FALSE,
    no_overlap = FALSE,
    cex = 1,
    cex_sub = 16 * cex,
    cex_point = 3 * cex,
    cex_lab = 19 * cex) {

    if (!isTRUE(text)) {
        func <- quote(geom_point(size = cex_point))
        if (!is.numeric(na.omit(group)))
            func$mapping <- aes(shape = as.factor(group))
    } else {

        f <- "geom_text"
        func <- quote(
            get(f)(aes(label = rownames(df)),
            size = cex_point)
        )
        
        if (no_overlap && NROW(df) <= 200) {
            f = paste0(f, '_repel')
            func$force = 0.2
            func$max.iter = 500
        }
    }

    if (is.null(p))
        p <- ggplot(df, aes(df[, 1], df[, 2], colour = as.factor(group)))

    if (length(name_group) > 15)
        name_group <- name_group[seq(15)]

    if (is.null(name_group))
        name_group <- 0
    
    axis <- function(margin){
        element_text(
            face = "italic",
            size = cex_lab * 0.75,
            margin = margin
        )
    }

    p <- p + eval(as.call(func)) + theme_classic() + geom_vline(
            xintercept = 0,
            col = "grey",
            linetype = "dashed",
            size = 1
        ) + geom_hline(
            yintercept = 0,
            col = "grey",
            linetype = "dashed",
            size = 1
        ) + labs(
                title = paste(title, "space"),
                x = print_comp(rgcca, compx, i_block),
                y = print_comp(rgcca, compy, i_block_y),
            color = name_group,
            shape = name_group
        ) + 
        scale_y_continuous(breaks = NULL) +
        scale_x_continuous(breaks = NULL) +
        theme_perso(cex, cex_sub) +
        theme(
            legend.key.width = unit(nchar(name_group), "mm"),
            axis.text = element_blank(),
            axis.title.y = axis(margin(0, 20, 0, 0)),
            axis.title.x = axis(margin(20, 0, 0, 0)),
            axis.line = element_blank()
        )

    if (length(unique(group)) != 1 && title == "Variable") {
        order_color(rgcca$a, p, collapse = collapse)
        # For qualitative response OR no response
    } else if ( is.character2(group[!is.na(group)]) ||
                length(unique(group)) <= 5 || 
            all( levels(as.factor(group)) %in% c("obs", "pred") )
        ) {
        p + scale_color_manual(values = color_group(group))
        # quantitative response
    } else
        p + scale_color_gradientn(colours = colours, na.value = "black")

}
