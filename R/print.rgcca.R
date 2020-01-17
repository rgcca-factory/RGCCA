print.rgcca <- function(x) 
{
  cat("Call: ")
  dput(x$call)
  cat("\n\n")
  
  if(is.list(x$crit))
  {
      critByNcomp=sapply(x$crit,function(t){return(t[length(t)])})
      cat("Sum_{j,k} c_jk g(cov(X_ja_j, X_ka_k) = ", sep = "", 
          paste(round(sum(critByNcomp), 4), sep = "", " "), fill = TRUE)
  }
  else
  {
      cat("Sum_{j,k} c_jk g(cov(X_ja_j, X_ka_k) = ", sep = "", 
          paste(round(x$crit[length(x$crit)], 4), sep = "", " "), fill = TRUE)
  }
  cat("There are J =", NCOL(x$C), "blocks.", fill = TRUE)
  cat("The design matrix is:\n") 
  colnames(x$C) = rownames(x$C) = names(x$a) ; print(x$C)
  cat("The", x$scheme, "scheme was used.", fill = TRUE)
  for (i in 1:NCOL(x$C)) {
    cat("The shrinkage parameter used for block", i, "was:", 
          round(x$tau[i], 4), fill = TRUE)
  }
}