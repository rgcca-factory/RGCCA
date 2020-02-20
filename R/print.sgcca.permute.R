#' Prints the results of permutation rgcca
#' @param x result of rgcca_permutation
print.permutation <- function (x) 
{
  cat("Call: ")
  dput(x$call)
  cat("\n\n")
  
  tab <- round(cbind(x$pvals, x$zstat, x$crit, rowMeans(x$permcrit)), 
               3)
  dimnames(tab) <- list(paste("Tuning parameter set ", 
                              sep = "", 1:length(x$pvals)), c("P-Value", 
                                                              "Z", "Crit", "Crit Perm"))
  print(tab, quote = FALSE)
  cat("Highest z score: ", max(x$zstat), "\n")
  cat("P-value corresponding to highest z score: ", x$pvals[which.max(x$zstat)], 
      fill = TRUE)
  cat("Tuning parameters corresponding to highest z score: ", 
      round(x$bestpenalties, 3), "\n")
  c1s <- round(x$penalties, 4)
  rownames(c1s) = 1:NROW(c1s)
    cat(fill = TRUE)
  cat("Tuning parameters used: ", fill = TRUE)
  print(c1s, quote = FALSE)
}