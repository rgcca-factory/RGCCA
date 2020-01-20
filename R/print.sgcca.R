print.sgcca <- function(x) 
{
  cat("Call: ")
  dput(x$call)
  cat("\n\n")
  cat("Sum_{j,k} c_jk g(cov(X_ja_j, X_ka_k) = ", sep = "", 
      paste(round(x$crit, 4), sep = "", " "), fill = TRUE)
  cat("There are J =", NCOL(x$call$C), "blocks.", fill = TRUE)
  cat("The design matrix is:\n") 
  colnames(x$call$C) = rownames(x$call$C) = names(x$a) ; print(x$call$C)
  cat("The", x$call$scheme, "scheme was used.", fill = TRUE)
  for (i in 1:NCOL(x$call$C)) {
    cat("\n The sparsity parameter used for block", i, "was:", 
          round(x$c1[i], 4), fill = TRUE)
    cat("Number of non-zero elements of canonical variate(s) for block ", 
        i, ": ", sep = "")
    if (is.matrix(x$a[[i]])) 
      cat(apply(x$a[[i]] != 0, 2, sum), fill = TRUE)
    if (!is.matrix(x$a[[i]])) 
      cat(sum(x$a[[i]] != 0), fill = TRUE)
    cat(fill = TRUE)
  }
}