library(plotly)

zstat <- as.data.frame(sgcca.res$zstat[, -1])
colnames(zstat) <- c(paste0("block_", 2:4), "z")

plotly::plot_ly(
    zstat,
    x = ~ block_2,
    y = ~ block_3,
    z = ~ block_4,
    marker = list(
        color = ~ z,
        colorscale = c('#FFE1A1', '#683531'),
        showscale = TRUE
    )
) %>% add_markers()
