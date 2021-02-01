#' Print a fitted rgcca_permutation object
#' 
#' Print a fitted rgcca_permutation object
#' @param x A fitted rgcca_permutation object (see  \code{\link[RGCCA]{rgcca_permutation}})
#' @param ... additional print parameters
#' @export
#' @examples
#' data(Russett)
#' A = list(agriculture = Russett[, seq(3)], 
#'          industry = Russett[, 4:5], 
#'          politic = Russett[, 6:11])
#'          
#' perm.out = rgcca_permutation(A, par_type = "tau", 
#'                              n_perms = 5, n_cores = 1)
#' print(perm.out)
print.permutation <- function(x, ...){
  
  cat("Call: ") 
  names_call <- c("type", "par_type", "n_perms", 
                  "quiet", "method", "tol", "scale", 
                  "scale_block", "superblock")
  char_to_print <- ""
  for (name in names_call) {
    if (name == "ncomp") {
      if (length(x$call$ncomp) > 1) {
        value <- (paste(x$call$ncomp, sep = "", collapse = ","))
        value <- paste0("c(", value, ")")
      }
    }
    if (name != "ncomp") {
      value <- x$call[[name]]
    }
    
    quo <- ifelse(is.character(value) & name != "ncomp", "'", "")
    vir <- ifelse(name == names_call[length(names_call)], "", ", ")
    char_to_print <- paste(char_to_print, name, "=", quo, value, 
                           quo, vir, collapse = "", sep = "")
  }
  cat(char_to_print, "\n")
  
  cat("The design matrix is:\n")
  colnames(x$call$connection) <- rownames(x$call$connection) <- names(x$a)
  print(x$call$connection)
  cat("\n")
  
  if (is.function(x$call$scheme)) {
    cat("The", deparse(x$call$scheme), "scheme was used.", fill = TRUE)
  } else {
    cat("The", x$call$scheme, "scheme was used.", fill = TRUE)
  }
  
  c1s <- round(x$penalties, 4)
  rownames(c1s) <- 1:NROW(c1s)
  cat(fill = TRUE)
  cat("Tuning parameters used: ", fill = TRUE)
  print(c1s, quote = FALSE, ...)
  cat("\n")
  
  tab <- round(cbind(x$crit, x$means, x$sds, x$zstat, x$pvals), 3)
  dimnames(tab) <- list(paste("Tuning parameter set ", sep = "", 
                              1:length(x$pvals)), 
                              c("crit", "crit perm", "sd", "zstat", "p-value"))
  print(tab, quote = FALSE, ...)
  
  cat(paste0("\nThe best combination was: ", 
             paste(round(x$bestpenalties, 3), collapse = ", "), 
             " for a z score of ", round(max(x$zstat), 3), 
             " and a p-value of ", round(x$pvals[which.max(x$zstat)], 3), 
             ".\n"))
}

