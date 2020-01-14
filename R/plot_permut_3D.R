#' Plot permuation in 3D
#' 
#' Plot permuation in 3D
#' 
#' @inheritParams plot3D
#' @inheritParams plot_permut_2D
#' @param sign A boolean to color by groups of alpha = 0.05, 0.01 or 0.001
#' @examples
#' data("Russett")
#' A = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
#'     politic = Russett[, 6:11] )
#' perm <- rgcca_permutation(A, nperm = 2, n_cores = 1)
#' plot_permut_3D(perm)
#' perm <- rgcca_permutation(A, p_c1 = TRUE, nperm = 2, n_cores = 1)
#' plot_permut_3D(perm)
# c1s <- expand.grid(
#     lapply(
#         seq(length(A)),
#         function(x) seq(1 / sqrt(ncol(A[[x]])), 1, by = 0.1)
#     )
# )
# perm <- rgcca_permutation(A, p_c1 = c1s, nperm = 2, n_cores = 1)
# plot_permut_3D(perm)
#' @export
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

    match.arg(type, c("crit", "zstat"))

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

    plotly::plot_ly(
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
            cauto = F,
            cmin = 0,
            cmax = max(zstat$z)
        )
    ) %>% 
    add_trace(type = "scatter3d", mode = "markers") %>% 
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
