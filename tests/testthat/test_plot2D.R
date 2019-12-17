
#' df = as.data.frame(matrix(runif(20*2, min = -1), 20, 2))
#' AVE = lapply(seq(4), function(x) runif(2))
#' rgcca_out = list(AVE = list(AVE_X = AVE))
#' plot2D(rgcca_out, df, "Samples", rep(c("a","b"), each=10), "Response")
