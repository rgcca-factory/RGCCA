color_nodes <- function(nodes) {
    unlist(lapply(as.list(1 - nodes$P / max(nodes$P)), 
        function(x)
            rgb(colorRamp( c("coral3", "khaki2"))(x) / 255)))
}
