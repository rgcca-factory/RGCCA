#'print.rgcca
#' Print rgcca results
#' @param object a result of rgcca function
#' @param ... other parameters used in print (for the displaying of matrices)
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

summary.rgcca <- function(object,...) 
{
  cat("Call: ")
  dput(object$call[!names(object$call)%in%c("blocks","connection")])
  cat("\n\n")
 
  if(is.list(object$crit))
  {
      critByNcomp=sapply(object$crit,function(t){return(t[length(t)])})
      cat("Sum_{j,k} c_jk g(cov(X_ja_j, X_ka_k) = ", sep = "", 
          paste(round(sum(critByNcomp), 4), sep = "", " "), fill = TRUE)
  }
  else
  {
      cat("Sum_{j,k} c_jk g(cov(X_ja_j, X_ka_k) = ", sep = "", 
          paste(round(object$crit[length(object$crit)], 4), sep = "", " "), fill = TRUE)
  }
  cat("There are J =", NCOL(object$call$connection), "blocks.", fill = TRUE)
  cat("The design matrix is:\n") 
  colnames(object$call$connection) = rownames(object$call$connection) = names(object$a) ; print(object$call$connection)
  if(is.function(object$call$scheme))
  {
      cat("The", deparse(object$call$scheme), "scheme was used.", fill = TRUE)
  }
  else
  {
      cat("The", object$call$scheme, "scheme was used.", fill = TRUE)
      
  }
  if(object$call$type %in% c("rgcca"))   
  {param="regularization"
    if(!is.matrix(object$tau))
    {
        for (i in 1:NCOL(object$call$connection)) {
            tau <- object$call$tau[i]
            if (tau != "optimal")
                tau <- round(tau , 4)
            
            

            cat("The",param," parameter used for block", i, "was:", 
                tau, fill = TRUE)
        }
    }
    else
    {
        cat("The",param," parameter used for block", i, "was: \n")
        print(round(object$tau,4),...)
    }
  }
  if(object$call$type %in% c("sgcca"))
      if(!is.matrix(object$sparsity))
      {
          for (i in 1:NCOL(object$call$connection)) {
              sparsity <- object$call$sparsity[i]
             param="shrinkage"
              cat("The",param," parameter used for block", i, "was:", 
                  sparsity, fill = TRUE)
          }
      }
  else
  {
      cat("The",param," parameter used for block", i, "was: \n")
      print(round(object$tau,4),...)
  }
}
