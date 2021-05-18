#' test_orderdf
 df = sapply(seq(2), function(x) runif(10))
 order_df(df)


 df = sapply(seq(2), function(x) runif(10))
 df[1, 2]=NA
 df[2,]=NA
 order_df(df)
