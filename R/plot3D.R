# Plot in 3 dimensions
# 
# Plot in 3 dimensions either to visualize the components of an analyse or the variables
# @inheritParams plot_ind
# @inheritParams plot2D
# @inheritParams get_comp
# @param type A character for the type of plot : either "ind" for individual plot or "var" for corcircle
plot3D <- function(
    df,
    rgcca_res,
    compx = 1,
    compy = 2,
    compz = 3,
    i_block = 1,
    i_block_y = i_block,
    i_block_z = i_block,
    text = TRUE,
    title = "Sample plot",
    colors = NULL,
    type = "ind",
    cex = 1,
    cex_point = 3 * cex,
    cex_lab = 19 * cex) {

    stopifnot(is(rgcca_res, "rgcca"))
    check_boolean("text", text)
    if (is.null(colors)) {
        colors <- "#A50026"
        colors[2] <- NA
        colors[3] <- "#313695"
    }else
        check_colors(colors)
    title <- paste0(title, collapse = " ")
    match.arg(type, c("ind", "var"))
    for (i in c("i_block", "i_block_y", "i_block_z"))
            check_blockx(i, get(i), rgcca_res$call$blocks)
    check_ncol(rgcca_res$Y, i_block)
    for (i in c("compx", "compy", "compz")) {
        if (!is.null(get(i)))
            check_compx(i, get(i), rgcca_res$call$ncomp, i_block)
    }
    for (i in c("cex", "cex_point", "cex_lab"))
        check_integer(i, get(i))

    load_libraries("plotly")
    `%>%` <- plotly::`%>%`

    if (is.na(colors[2]) && length(unique(df$resp)) == 1) {
        df$resp <- as.factor(rep("a", length(df$resp)))
        midcol <- "#cd5b45"
    } else
        midcol <- "gray"

    axis <- function(x, i)
        list(
                title = paste0("<i>", print_comp(rgcca_res, x, i), "</i>"),
                titlefont = list(
                        size = cex_lab * 0.75
                    )
            )

   if(!is.character(df$resp)){
      if (is.null(colors)){
        colors <- "#A50026"
        colors[2] <- NA
        colors[3] <- "#313695"
      }
   }

    color <- function(x){
        n <- length(x)
        if (!is.character(df$resp))
        { print(x)
          y=cut(
            x,
            breaks = n,
            labels = colorRampPalette(c(colors[1], colors[2], colors[3]))(n),
            include.lowest = TRUE)
          print(y)
          return(y)
        } 
        else
        {
          y=color_group(seq(length(unique(df$resp))), colors = colors)
          return(y)
        }
    }

    subdf <- function(x) 
        df[which(df$resp == levels(df$resp)[x]), ]

    add_trace_manual <- function(p, x){

        l <- levels(df$resp)

        func <- quote(
            plotly::add_trace(
                p,
                name = l[x],
                x = ~ subdf(x)[, 1],
                y = ~ subdf(x)[, 2],
                z = ~ subdf(x)[, 3],
                type = "scatter3d",
                showlegend = TRUE
            )
        )

        color <- color_group(seq(length(l)), colors = colors)[x]

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

    if (!is.character(df$resp)) {

        if (text)
            visible <- "legendonly"
        else
            visible <- TRUE

        p <- plotly::plot_ly(
            name = "samples",
            x = ~ df[, 1],
            y = ~ df[, 2],
            z = ~ df[, 3],
            mode = "markers",
            type = "scatter3d",
            showlegend = FALSE,
            color = df$resp,
            size = I(200),
            colors = c(colors[1], colors[2], colors[3]),
            visible = visible
        )

        if (text) {
          print("in")
            p <- p %>%
              plotly::add_trace(
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
        p <- plotly::plot_ly()
        
        for (i in seq(length(levels(df$resp))))
            p <- p %>% add_trace_manual(i)
    }

    p <- p %>%
      plotly::layout(
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
        plotly::add_trace(
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
