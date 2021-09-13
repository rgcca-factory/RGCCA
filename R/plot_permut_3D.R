# Plot permuation in 3D
#
# Plot permuation in 3D
#
# @inheritParams plot3D
# @inheritParams plot_permut_2D
# @param sign A boolean to color by groups of alpha = 0.05, 0.01 or 0.001
#@export
plot_permut_3D <- function(
    perm,
    type = "zstat",
    sign = FALSE,
    i_block = 1,
    i_block_y = 2,
    i_block_z = 3,
    cex = 1,
    cex_point = 3 * cex,
    cex_lab = 19 * cex) {

    stopifnot(is(perm, "permutation"))
    match.arg(type, c("crit", "zstat"))
    check_boolean("sign", sign)
    for (i in c("cex_point", "cex_lab"))
        check_integer(i, get(i))
    check_integer("cex", cex, float = TRUE)
    for (i in c("i_block", "i_block_y", "i_block_z"))
        check_blockx(i, get(i), perm$penalties[1,])

    load_libraries("plotly")
    `%>%` <- plotly::`%>%`

    switch(
        type,
        "zstat" =  y_title <- "Z-score",
        "crit" = y_title <- "RGCCA criterion"
    )

    zstat <- as.data.frame(
        cbind(perm$penalties[,c(i_block, i_block_y, i_block_z)],
        z = unlist(perm[type])))
    best <- which.max(zstat$z)

    axis3D_perm <- function(i)
        axis3D(colnames(zstat)[i], cex_lab)

    if (type == "zstat" && sign)
        zstat[, 4] <- as.double(zstat$z > qnorm(1 - 0.05 / 2))
            + as.double(zstat$z > qnorm(1 - 0.01 / 2))
            + as.double(zstat$z > qnorm(1 - 0.001 / 2))

    p <- plotly::plot_ly(
        zstat,
        x = ~ zstat[, 1],
        y = ~ zstat[, 2],
        z = ~ zstat[, 3],
        marker = list(
            color = ~ zstat[, 4],
            size = cex_point * 2,
            showscale = TRUE,
            colorbar = list(title = y_title),
            colorscale = list(
                list(0, "rgb(165,0,38)"),
                list(mean(zstat$z), "rgb(254,224,144)"),
                list(max(zstat$z), "rgb(49,54,149)")
            ),
            cauto = FALSE,
            cmin = 0,
            cmax = max(zstat$z)
        )
    )
    plotly::add_trace(p, type = "scatter3d", mode = "markers") %>%
    layout3D(
        title = paste0(
            "Permutation scores \n(best value : ",
            paste(round(perm$penalties[best,], 2), collapse = ", "),
            ")"
        ),
        cex,
        scene = list(
            xaxis = axis3D_perm(1),
            yaxis = axis3D_perm(2),
            zaxis = axis3D_perm(3)
        ))

}
