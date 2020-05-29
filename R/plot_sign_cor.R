plot_sign_cor <-  function(df, compx = 1, n_mark = 35){
  if(n_mark < NROW(df))
    df = tail(df, n_mark)
  df = data.frame(x = log10(df[, compx]), order = seq(NROW(df)), row.names = row.names(df), color = log10(df[, compx]))
  p = ggplot(df, aes(order, x, fill = color))
  plot_histogram(p, df, "", 1, cex_sub = 25, cex_axis = 15)  +
    labs(fill = "Log10 \np-adjusted")
}