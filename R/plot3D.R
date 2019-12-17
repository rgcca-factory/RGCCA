#' Plot in 3 dimensions
#' 
#' Plot in 3 dimensions either to visualize the components of an analyse or the variables
#' @inheritParams plot_ind
#' @inheritParams plot2D
#' @inheritParams get_comp
#' @param type A character for the type of plot : either "ind" for individual plot or "var" for corcircle
#' @examples
#' library(RGCCA)
#' data("Russett")
#' blocks = list(agriculture = Russett[, seq(3)],
#'     politic = Russett[, 6:11] )
#' rgcca_out = rgcca.analyze(blocks, ncomp = rep(3, 2))
#' df = get_comp(rgcca_out, compz = 3)
#' plot3D(df, rgcca_out, i_block = 2)
#' plot3D(df, rgcca_out, i_block = 2, text = FALSE)
#' response = factor( apply(Russett[, 9:11], 1, which.max),
#'                   labels = colnames(Russett)[9:11] )
#' response = blocks[[2]][, 1]
#' names(response) = row.names(blocks[[2]])
#' df = get_comp(rgcca_out, response, compz = 3)
#' plot3D(df, rgcca_out, i_block = 2, text = FALSE)
#' plot3D(df, rgcca_out, i_block = 2)
#' df = get_ctr2(rgcca_out, compz = 3, i_block = 1, collapse = TRUE)
#' plot3D(df, rgcca_out, i_block = 2, type = "var")
#' @export
plot3D <- function(
    df,
    rgcca,
    compx = 1,
    compy = 2,
    compz = 3,
    i_block = 1,
    i_block_y = i_block,
    i_block_z = i_block,
    text = TRUE,
    title = "Sample plot",
    type = "ind",
    cex = 1,
    cex_point = 3 * cex,
    cex_lab = 19 * cex) {

    if (length(unique(df$resp)) == 1) {
        df$resp = as.factor(rep("a", length(df$resp)))
        midcol = "#cd5b45"
    } else
        midcol = "gray"

    axis <- function(x, i)
        list(
                title = paste0("<i>", print_comp(rgcca, x, i), "</i>"),
                titlefont = list(
                        size = cex_lab * 0.75
                    )
            )

    color <- function(x){
        n <- length(x)
        if (!is.character2(df$resp))
            cut(
                x,
                breaks = n,
                labels = colorRampPalette(c("#A50026", midcol,  "#313695"))(n),
                include.lowest = TRUE)
        else
            color_group(seq(length(unique(df$resp))))
    }

    subdf <- function(x) 
        df[which(df$resp == levels(df$resp)[x]), ]

    add_trace_manual <- function(p, x){

        l <- levels(df$resp)

        func <- quote(
            add_trace(
                p,
                name = l[x],
                x = ~ subdf(x)[, 1],
                y = ~ subdf(x)[, 2],
                z = ~ subdf(x)[, 3],
                type = "scatter3d",
                showlegend = TRUE
            )
        )

        color <- color_group(seq(length(l)))[x]

        if (text) {
            func$mode <- "text"
            func$text <- ~row.names(subdf(x))
            func$textfont <- list(
                color = color,
                size = cex_point * 4
            )
        }else{
            func$mode <- "markers"
            func$marker <- list(
                color = color,
                size = cex_point * 1.5
            )
        }

        eval(func)
    }


    if (!is.character2(df$resp)) {

        if (text)
            visible <- "legendonly"
        else
            visible <- TRUE

        p <- plot_ly(
            name = "samples",
            x = ~ df[, 1],
            y = ~ df[, 2],
            z = ~ df[, 3],
            mode = "markers",
            type = "scatter3d",
            showlegend = FALSE,
            color = df$resp,
            size = I(200),
            colors = c("#A50026", midcol,  "#313695"),
            visible = visible
        )

        if (text) {
            p <- p %>%
                add_trace(
                    name = "samples",
                    x = ~ df[, 1],
                    y = ~ df[, 2],
                    z = ~ df[, 3],
                    mode = "text",
                    type = "scatter3d",
                    text = ~ row.names(df),
                    textfont = list(
                        color = color(df$resp),
                        size = cex_point * 4
                    ),
                    showlegend = FALSE,
                    visible = TRUE
                )
        }

    }else{
        p <- plot_ly()
        
        for (i in seq(length(levels(df$resp))))
            p <- p %>% add_trace_manual(i)
    }

    p <- p %>%
        layout(
            autosize = TRUE,
            margin = list(
                l = 50,
                r = 50,
                b = 50,
                t = 100
            ),
            scene = list(
                aspectmode = 'cube',
                xaxis = axis(compx, i_block),
                yaxis = axis(compy, i_block_y),
                zaxis = axis(compz, i_block_z)
            ),
            title = list(
                text = paste0('<b>', title, '</b>'),
                font = list(
                    size = 25 * cex,
                    face = "bold"
                )
            )
        )

    plot_circle3D <- function(p, x, y, z){
        df <- cbind(plot_circle(), 0)
        add_trace(
            p = p,
            x = df[, x],
            y = df[, y],
            z = df[, z],
            showlegend = FALSE,
            hoverinfo = "none" ,
            mode = "lines",
            type = "scatter3d",
            line = list(color = "grey", width = 4)
        )
    }

    if (type == "var")
        p <- p %>% plot_circle3D(1, 2, 3) %>% plot_circle3D(1, 3, 2) # %>% plot_circle3D(3, 2, 1)

    return(p)
}
