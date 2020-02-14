#'print.rgcca
#' Print rgcca results
#' @param x a result of rgcca function
#' @param ... other parameters 
#' @export
#' @examples 
#' data(Russett)
#'X_agric =as.matrix(Russett[,c("gini","farm","rent")]);
#'X_ind = as.matrix(Russett[,c("gnpr","labo")]);
#'X_polit = as.matrix(Russett[ , c("demostab", "dictator")]);
#'A = list(X_agric, X_ind, X_polit);
#'C = matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3);
#'res = rgcca(A, C, ncomp=rep(2,3),tau = c(1, 1, 1), scheme = "factorial", scale = TRUE,verbose=FALSE)
#'print(res)

print.rgcca <- function(x,...) 
{
  cat("Call: ")
  dput(x$call[!names(x$call)%in%c("blocks","connection")])
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
  cat("There are J =", NCOL(x$call$connection), "blocks.", fill = TRUE)
  cat("The design matrix is:\n") 
  colnames(x$call$connection) = rownames(x$call$connection) = names(x$a) ; print(x$call$connection)
  if(is.function(x$call$scheme))
  {
      cat("The", deparse(x$call$scheme), "scheme was used.", fill = TRUE)
  }
  else
  {
      cat("The", x$call$scheme, "scheme was used.", fill = TRUE)
      
  }
  if(x$call$type %in% c("rgcca","sgcca"))   
  {
    if(!is.matrix(x$tau))
    {
        for (i in 1:NCOL(x$call$connection)) {
            tau <- x$call$tau[i]
            if (tau != "optimal")
                tau <- round(tau , 4)
            if(x$call$type=="rgcca"){param="regularization"}
            if(x$call$type=="sgcca"){param="shrinkage"}
            cat("The",param," parameter used for block", i, "was:", 
                tau, fill = TRUE,...)
        }
    }
    else
    {
        cat("The",param," parameter used for block", i, "was: \n")
        print(round(x$tau,4))
    }
    
  }
}
