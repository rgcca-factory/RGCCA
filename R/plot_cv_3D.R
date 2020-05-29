#' Plot cross-validation in 3D
#' 
#' Plot cross-validation in 3D
#' 
#' @inheritParams plot3D
#' @inheritParams plot_permut_2D
# @export
plot_cv_3D <- function(
  scores_res,
  i_block = 1,
  i_block_y = 2,
  i_block_z = 3,
  cex = 1,
  cex_point = 3 * cex,
  cex_lab = 19 * cex,
  mean_col = quantile(scores_res$scores, probs = 0.25),
  max_col = quantile(scores_res$scores, probs = 0.5),
  min_col = min(scores_res$scores)) {
  
  # stopifnot(is(perm, "permutation"))
  for (i in c("cex","cex_point", "cex_lab"))
    check_integer(i, get(i))
  for (i in c("i_block", "i_block_y", "i_block_z"))
    check_blockx(i, get(i), scores_res[1, -1])
  
  load_libraries("plotly")
  `%>%` <- plotly::`%>%`
  
  axis3D_perm <- function(i)
    axis3D(colnames(scores_res)[i], cex_lab)
  
  p <- plotly::plot_ly(
    scores_res,
    x = ~ scores_res[, 1],
    y = ~ scores_res[, 2],
    z = ~ scores_res[, 3],
    marker = list(
      color = ~ scores_res[, 4],
      size = cex_point * 2,
      showscale = TRUE,
      colorbar = list(title = "RMSE"),
      colorscale = list(
        list(max_col, "rgb(49,54,149)"), 
        list(mean_col, "rgb(254,224,144)"),
        list(min_col, "rgb(165,0,38)")
      ),
      cauto = F,
      cmin = max_col,
      cmax = min_col
    )
  )
  plotly::add_trace(p, type = "scatter3d", mode = "markers") %>% 
    layout3D(
      title = paste0(
        "Cross-validation \n(best value : ",
        paste(round(scores_res[which.min(scores_res$scores), -4], 2), collapse = ", "),
        ")"
      ),
      cex,
      scene = list(
        xaxis = axis3D_perm(1),
        yaxis = axis3D_perm(2),
        zaxis = axis3D_perm(3)
      ))
  
}