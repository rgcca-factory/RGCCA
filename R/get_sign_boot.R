get_sign_boot <- function(df_b, x = 0.05) {
  n_boot <- attributes(df_b)$n_boot
  nvar <- attributes(df_b)$n_var
  avg_n_occ <- sum(df_b$occurrences*n_boot) / n_boot
  p <- qbinom(size = n_boot, prob = avg_n_occ / nvar, p = 1 - x / nvar) / n_boot
  df_b <- order_df(df_b[, -NCOL(df_b)], "occurrences")
  row.names(df_b[df_b$occurrences > p, ])
}