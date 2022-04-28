#' Print rgcca
#
#' @title Print the call of rgcca results
#' @param x A RGCCA object (see \code{\link{rgcca}})
#' @param ... other parameters used in print (for the displaying of matrices)
#' @export
#' @examples
#' data(Russett)
#' X_agric <- as.matrix(Russett[, c("gini", "farm", "rent")])
#' X_ind <- as.matrix(Russett[, c("gnpr", "labo")])
#' X_polit <- as.matrix(Russett[, c("demostab", "dictator")])
#' A <- list(X_agric, X_ind, X_polit)
#' C <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3)
#' res <- rgcca(A,
#'   connection = C, ncomp = rep(2, 3), tau = c(1, 1, 1),
#'   scheme = "factorial", scale = TRUE, verbose = FALSE
#' )
#' print(res)
print.rgcca <- function(x, ...) {
  ### Print parameters of the function
  print_call(x$call)

  ### Print criterion
  if (is.list(x$crit)) {
    crit <- Reduce("+", lapply(x$crit, function(t) {
      return(t[length(t)])
    }))
  } else {
    crit <- x$crit[length(x$crit)]
  }
  cat("Sum_{j,k} c_jk g(cov(X_ja_j, X_ka_k) = ",
    sep = "",
    paste(round(crit, 4), sep = "", " "), fill = TRUE
  )

  ### Print regularization parameter or the number of selected variables
  cat("\n")
  if (!tolower(x$call$method) %in% c("sgcca", "spca", "spls")) {
    param <- "regularization"
    if (!is.matrix(x$call$tau)) {
      for (i in seq(NCOL(x$call$connection))) {
        tau <- x$call$tau[i]
        cat("The", param, "parameter used for", names(x$call$blocks)[i],
          "is:", round(tau, 4),
          fill = TRUE
        )
      }
    } else {
      cat("The", param, "parameters used are: \n")
      print(round(x$call$tau, 4), ...)
    }
  } else {
    nb_selected_var <- lapply(
      x$a,
      function(a) apply(a, 2, function(l) sum(l != 0))
    )
    param <- "sparsity"
    if (!is.matrix(x$call$sparsity)) {
      for (i in seq(NCOL(x$call$connection))) {
        sparsity <- x$call$sparsity[i]

        cat("The", param, "parameter used for", names(x$call$blocks)[i], "is:",
          sparsity, "(with", paste(nb_selected_var[[i]], collapse = ", "),
          "variables selected)",
          fill = TRUE
        )
      }
    } else {
      cat("The", param, "parameters used are: \n")
      print(round(x$call$sparsity, 4), ...)
      cat("The number of selected variables are: \n")
      print(do.call(cbind, nb_selected_var))
    }
  }
}
