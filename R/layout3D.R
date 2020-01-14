layout3D <- function(p, title, cex, scene){
    layout(
        p,
        autosize = TRUE,
        margin = list(
            l = 50,
            r = 50,
            b = 50,
            t = 100
        ),
        scene = scene,
        title = list(
            text = paste0('<b>', title, '</b>'),
            font = list(size = 25 * cex)
        )
    )
}
