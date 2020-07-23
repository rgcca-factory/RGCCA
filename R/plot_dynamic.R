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
plot_dynamic <- function(
    f,
    ax = NULL,
    text = "name+x+y",
    dynamicTicks = FALSE,
    type = "regular") {

    if (is.null(ax))
        ax <- list(linecolor = "white",
            ticks = "",
            titlefont = list(size = 23))

    `%>%` <- plotly::`%>%`
    # Convert a ggplot into a plotly object add a layout with predefined
    # formats for
    # x- and y- axis set the style to show onMouseOver text
    if (type == "regular")
        p <- plotly::plotly_build(
            plotly::ggplotly(f, dynamicTicks = dynamicTicks) %>%
                plotly::layout(
                    xaxis = ax,
                    yaxis = ax,
                    annotations = list(showarrow = FALSE, text = "")
                ) %>% 
                plotly::style(hoverinfo = text)
            )
    else 
        p <- ggplotly(f)

    legend_qual <- p$x$layout$annotations[[1]]$text
    legend_quant <- p$x$data[[length(p$x$data)]]$marker$colorbar$title
    if (!is.null(legend_qual))
        legend_title <- legend_qual
    else
        legend_title <- legend_quant

    if (!is.null(legend_title) && legend_title != "") {
        # Deals with a too short name of modalities
        p$x$layout$margin$r <- max(nchar(strsplit(legend_title, "\n")[[1]])) * 13
        # if \n in title
        if (type != "boot1D")
            legend_title <- strsplit(legend_title, "\n")[[1]]

        # for shiny corcircle, if text = TRUE, two legends will appear.
        # only the first one will be selected
        legend_title <- unlist(strsplit(legend_title, "<br />"))[1]

        # to prevent print a 'NA' when there is no legend in plot
        if (is.na(legend_title))
            legend_title <- ""

        # set the font for this title
        if (!is.null(legend_qual)) {
            p$x$layout$annotations[[1]]$text <- paste0("<i><b>", legend_title, "</b></i>")
            # set on the top the position of the legend title
            p$x$layout$annotations[[1]]$yanchor <- "top"
                if (length(legend_title) > 1)
            p$x$layout$annotations[[1]]$y = 1.05
        } else
            p$x$data[[length(p$x$data)]]$marker$colorbar$title <- paste0("<i><b>", legend_title, "</b></i>")
    }

    if (!is.null(f$labels$subtitle)) {
        if (packageVersion("plotly") < 4.9)
            p$x$layout$title <- paste0(
                    p$x$layout$title,
                    "<br><i>",
                    "c",
                    substring(f$labels$subtitle, 2),
                    "</i>"
                )
        else
            p$x$layout$title$text <- paste0(
                p$x$layout$title$text,
                "<br><i>",
                "c",
                substring(f$labels$subtitle, 2),
                "</i>"
            )
    }

    if (NCOL(f$data) == 3)
        p$sample_names <- lapply(
            levels(as.factor(f$data[, 3])), 
            function(x) row.names(subset(f$data, f$data[, 3] == x)))
    else
        p$sample_names <- list(row.names(f$data))

    p$x$layout$margin$t <- 100

    if (!grepl("1D", type) && !type %in% c("perm", "cv")) {
        p <- config(
            p,
            modeBarButtons = list(
                list("toImage"),
                list("resetScale2d"),
                list("zoom2d"),
                list("pan2d")
            ),
            displaylogo = FALSE
        )
    } else
        p <- config(p, modeBarButtons = list(list("toImage")))

    p$x$layout$title$y <- 0.95

    config(
        p,
        editable = TRUE,
        displaylogo = FALSE,
        edits = list(shapePosition = F)
    )  %>% 
        layout(hovermode = "closest")
}
