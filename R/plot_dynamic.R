# Dynamic visualization of the outputs
# f: ggplot2 function
# ax: list object containing attributes of xaxis / yaxis parameter in plotly
# (https://plot.ly/javascript/reference/, xaxis/yaxis)
# text: axis information to print (among y, x, x+y, text)
#  (https://plot.ly/javascript/reference/, hoverinfo)
# dynamicTicks: a bolean giving the generation of axis tick labels
# (otherwhise samplesPlot which do not have traces could not be convereted
# in ggplotly)
# return a plotly object
plot_dynamic <- function(f,
    ax = NULL,
    text = "name+x+y",
    legend = TRUE,
    dynamicTicks = FALSE) {
    
    if (is.null(ax))
        ax <- list(linecolor = "white",
            ticks = "",
            titlefont = list(size = 23))

    # Convert a ggplot into a plotly object add a layout with predefined
    # formats for
    # x- and y- axis set the style to show onMouseOver text
    p <- plotly_build(
            ggplotly(f, dynamicTicks = dynamicTicks) %>%
            layout(
                xaxis = ax,
                yaxis = ax,
                annotations = list(showarrow = FALSE, text = "")
            ) %>% 
            style(hoverinfo = text)
        )
    
    if (legend) {
        
        # set on the top the position of the legend title
        p$x$layout$annotations[[1]]$yanchor <- "top"
        # Deals with a too short name of modalities
        p$x$layout$margin$r <- nchar(p$x$layout$annotations[[1]]$text) * 13
        p$x$layout$margin$t <- 100
        # for shiny corcircle, if text = TRUE, two legends will appear.
        # only the first one will be selected
        title <- unlist(strsplit(p$x$layout$annotations[[1]]$text, "<br />"))[1]

        # to prevent print a 'NA' when there is no legend in plot
        if (is.na(title))
            title <- ""

        # set the font for this title
        p$x$layout$annotations[[1]]$text <- paste0("<i><b>", title, "</b></i>")
        # Sys.info()[['sysname']]

        if (!is.null(f$labels$subtitle)) {
            if (packageVersion("plotly") < 4.9)
                p$x$layout$title <- paste0(
                        p$x$layout$title,
                        "<br><b>",
                        "c",
                        substring(f$labels$subtitle, 2),
                        "</b>"
                    )
            else
                p$x$layout$title$text <- paste0(
                    p$x$layout$title$text,
                    "<br><b>",
                    "c",
                    substring(f$labels$subtitle, 2),
                    "</b>"
                )
        }
    }
    
    if (NCOL(f$data) == 3)
        p$sample_names <- lapply(
            levels(as.factor(f$data[, 3])), 
            function(x) row.names(subset(f$data, f$data[, 3] == x)))
    else
        p$sample_names <- list(row.names(f$data))
    return(p)
}
