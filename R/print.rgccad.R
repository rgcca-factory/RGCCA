print.rgccad <- function(x, ...) 
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
  cat("There are J =", NCOL(x$call$C), "blocks.", fill = TRUE)
  cat("The design matrix is:\n") 
  colnames(x$call$C) = rownames(x$call$C) = names(x$a) ; print(x$call$C)
  cat("The", x$call$scheme, "scheme was used.", fill = TRUE)
  
  x$call$tau <- elongate_arg(x$call$tau, x$a)
  for (i in 1:NCOL(x$call$C)) {
    tau <- x$call$tau[i]
      if (tau != "optimal")
          tau <- round(tau , 4)
      cat("The shrinkage parameter used for block", i, "was:", 
          tau, fill = TRUE,...)
  }
}
