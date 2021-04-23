#'# get_block_var
#'''
 df = data.frame(x = runif(30), order = 30:1)
 library("ggplot2")
 p = ggplot(df, aes(order, x))
 plot_histogram(p, df, "This is my title")
 # Add colors per levels of a variable
 df$color = rep(c(1,2,3), each=10)
 p = ggplot(df, aes(order, x, fill = as.character(color)))
 plot_histogram(p, df, "Histogram", as.character(df$color))
