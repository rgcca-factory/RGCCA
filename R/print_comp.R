#' Print the variance of a component
#'
#' Print the percent of explained variance for a component of a block
#' (by default, the superblock or the last one) analysed by R/SGCCA
#'
#' @inheritParams plot_ind
#' @param n An integer giving the index of the analysis component
#' @param i An integer giving the index of a list of blocks
#' @param outer A boolean for ave plot case
#' @return A character for the variance on the component
#' @seealso \code{\link[RGCCA]{rgccad}}, \code{\link[RGCCA]{sgcca}}

print_comp <- function(rgcca_res, n = 1, i = length(rgcca_res$AVE$AVE_X),
                       outer = FALSE) {
  nvar <- sum(rgcca_res$a[[i]][, n] != 0)
  if (
    !tolower(rgcca_res$call$method) %in% c("spls", "spca", "sgcca") ||
      nvar == length(rgcca_res$a[[i]][, n])
  ) {
    var_text <- ""
  } else {
    var_text <- paste0(nvar, " variables, ")
  }

  ave <- quote(paste0(round(AVE[n] * 100, 1), "%"))
  if (isTRUE(outer)) {
    AVE <- rgcca_res$AVE$AVE_outer
    if (length(rgcca_res$AVE$AVE_outer) > 1) {
      n <- seq(2)
    } else {
      n <- 1
    }
    paste0("First outer comp. : ", paste(eval(ave), collapse = " & "))
  } else {
    AVE <- rgcca_res$AVE$AVE_X[[i]]
    paste0("Comp. ", n, " (", var_text, eval(ave), ")")
  }
}
