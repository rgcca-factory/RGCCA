#' Histogram of Average Variance Explained
#'

 random_val = function(y=1) lapply(seq(4),
function(x) matrix(runif(4), y, 2))
rgcca_out = list(AVE = list(AVE_X = random_val()),
     a = random_val(2), ncomp = rep(2, 4))
names(rgcca_out$a) <- LETTERS[seq(4)]
library("ggplot2")
 for(i in seq(1,4))
 names(rgcca_out$AVE$AVE_X[[i]]) <- c(1,2)
plot_ave(rgcca_out)
