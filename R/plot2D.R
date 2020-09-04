# 'Plot of components space
#'
#' Plots RGCCA components in a bi-dimensional space
#'
#' @param compx Integer corresponding to the number of x-component
#' @param compy Integer corresponding to the number of y-component
#' @param i_block Integer corresponding to the first block
#' @param text If TRUE, labels are plotted. If FALSE points are plotted.
#' @param i_block_y Integer corresponding to the second block
#' @param no_overlap If TRUE, the potential overlaps are reduced
#' @param rgcca_res Result of rgcca function
#' @param df A dataframe
#' @param title A character with the name of the space (either "Variables" or
#' "Samples")
#' @param group A vector of character with levels used to color the points
#' @param name_group A character giving the type of groups (either "Blocs" or
# "Response")
#' @param p A ggplot object
#' @param colors A vectof of character to color quantitative data
#' @param cex An integer for the size of the plot parameters
#' @param cex_main An integer for the size of the title
#' @param cex_sub An integer for the size of the subtitle
#' @param cex_point An integer for the size of the points or the text in the plot
#' @param cex_lab An integer for the size of the axis titles
#' @param collapse A boolean to combine the variables of each block as result
#' @importFrom ggplot2 ggplot
# @examples
# df = as.data.frame(matrix(runif(20*2, min = -1), 20, 2))
# AVE = lapply(seq(4), function(x) runif(2))
# rgcca_out = list(AVE = list(AVE_X = AVE), call = list(type = "rgcca"))
# plot2D(rgcca_out, df, "Samples", rep(c("a","b"), each=10), "Response")
# data(Russett)
# blocks = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
# politic = Russett[, 6:11] )
# rgcca_out = rgcca(blocks)
# plot2D(rgcca_out, df)

plot2D <- function(
    rgcca_res,
    df,
    title = "",
    group = 1,
    name_group = "Response",
    compx = 1,
    compy = 2,
    i_block = 1,
    p = NULL,
    text = TRUE,
    i_block_y = i_block,
    colors = NULL,
    collapse = FALSE,
    no_overlap = FALSE,
    cex = 1,
    cex_main = 14 * cex,
    cex_sub = 12 * cex,
    cex_point = 3 * cex,
    cex_lab = 10 * cex) {
    
    title <- paste0(title, collapse = " ")
    name_group <- paste0(name_group, collapse = " ")
    check_colors(colors)
    for (i in c("cex", "cex_main", "cex_sub", "cex_point", "cex_lab"))
        check_integer(i, get(i))
    for (i in c("text", "no_overlap"))
        check_boolean(i, get(i))
    if (NROW(df) > 100)
        cex_point <- 2

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
            tryCatch(load_libraries("ggrepel"), silent = TRUE)
            if (("ggrepel" %in% installed.packages()[, "Package"])) {
                f <- paste0(f, '_repel')
                func$force <- 0.2
                func$max.iter <- 500
            }
            else{
                warning("Please install ggrepel.")
            }
        }
    }

    if (is.null(p)) {
        p <- ggplot(df, aes(df[, 1], df[, 2], colour = as.factor(group)))
    }

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
    names_blocks=names(rgcca_res$call$blocks)
    if(is.null(names_blocks))
    {
        names_blocks=paste("Block", 1:length(rgcca_res$call$blocks))
    }
    if(i_block!=i_block_y)
    {
        xlab=paste0(print_comp(rgcca_res, compx, i_block)," - ", names_blocks[i_block])
        ylab=paste0(print_comp(rgcca_res, compy, i_block_y)," - ", names_blocks[i_block_y])
    }
    else
    {
        xlab=print_comp(rgcca_res, compx, i_block)
        ylab=print_comp(rgcca_res, compy, i_block_y)
        
    }
        
    p <- p + eval(as.call(func)) + theme_classic() + geom_vline(
            xintercept = 0,
            col = "grey",
            linetype = "dashed",
            size = 0.5
        ) + geom_hline(
            yintercept = 0,
            col = "grey",
            linetype = "dashed",
            size = 0.5
        ) + labs(
                title = title,
                x = xlab,
                y = ylab,
            color = name_group,
            shape = name_group
        ) + 
        scale_y_continuous(breaks = NULL) +
        scale_x_continuous(breaks = NULL) +
        theme_perso(cex, cex_main, cex_sub) +
        theme(
            legend.key.width = unit(nchar(name_group), "mm"),
            axis.text = element_blank(),
            axis.title.y = axis(margin(0, 20, 0, 0)),
            axis.title.x = axis(margin(20, 0, 0, 0)),
            axis.line = element_blank()
        )

    
    if (length(unique(group)) != 1 && is(df, "d_var2D")) {
        p <- order_color(rgcca_res$a, p, collapse = collapse, colors = colors)
        # For qualitative response OR no response
    } else if ( is.character2(group[!is.na(group)]) ||
                length(unique(group)) <= 5 || 
            all( levels(as.factor(group)) %in% c("obs", "pred") )
        ) {
        p <- p + scale_color_manual(values = color_group(group, colors))
        # quantitative response
    } else{
        if (is.null(colors))
            colors <- c("blue", "gray", "#cd5b45")
        p <- p + scale_color_gradientn(colours = colors, na.value = "black")
    }

    return(p)
}
