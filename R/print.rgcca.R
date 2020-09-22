#' Print rgcca
#
#' @title Print the call of rgcca results
#' @param x A RGCCA object (see \code{\link{rgcca}})
#' @param ... other parameters used in print (for the displaying of matrices)
#' @export
#' @examples 
#' data(Russett)
#'X_agric =as.matrix(Russett[,c("gini","farm","rent")]);
#'X_ind = as.matrix(Russett[,c("gnpr","labo")]);
#'X_polit = as.matrix(Russett[ , c("demostab", "dictator")]);
#'A = list(X_agric, X_ind, X_polit);
#'C = matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3);
#'res = rgcca(A, connection=C, ncomp=rep(2,3),tau = c(1, 1, 1), 
#'scheme = "factorial", scale = TRUE,verbose=FALSE)
#'print(res)

print.rgcca <- function(x,...) 
{
  cat("Call: ")
  names_call=c("type","superblock","scale","scale_block","init","bias","tol","method","ncomp")
  char_to_print=""
  for(name in names_call)
  {
      if(name=="ncomp"){if(length(x$call$ncomp)>1){value=(paste(x$call$ncomp,sep="",collapse=","));value=paste0("c(",value,")")}}
      if(name!="ncomp"){value=x$call[[name]]}
      quo=ifelse(is.character(value)&name!="ncomp","'","")
      vir=ifelse(name==names_call[length(names_call)],"",", ")
      char_to_print=paste(char_to_print,name,'=',quo,value,quo,vir, collapse="",sep="")
  }
  cat(char_to_print)
   
  cat("\n")
  cat("There are J =", NCOL(x$call$connection), "blocks.", fill = TRUE)
  cat("The design matrix is:\n") 
  colnames(x$call$connection) = rownames(x$call$connection) = names(x$a) ; print(x$call$connection)
  cat("\n")
   if(is.function(x$call$scheme))
  {
      cat("The", deparse(x$call$scheme), "scheme was used.", fill = TRUE)
  }
  else
  {
      cat("The", x$call$scheme, "scheme was used.", fill = TRUE)
      
  }
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
  cat("\n")
  if(x$call$type %in% c("rgcca"))   
  {
     param="regularization"
    if(!is.matrix(x$tau))
    {
        for (i in 1:NCOL(x$call$connection)) 
        {
            tau <- x$call$tau[i]
             cat("The",param,"parameter used for block", i, "was:",  round(tau,4), fill = TRUE)
        }
    }
    else
    {
        
        cat("The",param,"parameters used were: \n")
        print(round(x$tau,4),...)
    }
  }
  if(x$call$type %in% c("sgcca"))
  {
      param="sparsity"
      if(!is.matrix(x$sparsity))
      {
          for (i in 1:NCOL(x$call$connection)) {
              sparsity <- x$call$sparsity[i]
              
              cat("The",param,"parameter used for block", i, "was:", 
                  sparsity, fill = TRUE)
          }
      }
      else
      {
          for (i in 1:NCOL(x$call$connection)) 
        {
          cat("The",param,"parameters used were: \n")
          print(round(x$sparsity,4),...)
         }
      }
  }
   
}
