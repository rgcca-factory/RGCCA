#' Print the variance of a component
#'
#' Prints the percent of explained variance for a component of a block 
#' (by default, the superblock or the last one) analysed by R/SGCCA
#'
#' @inheritParams plot_ind
#' @param n An integer giving the index of the analysis component
#' @param i An integer giving the index of a list of blocks
#' @param outer A boolean for ave plot case
#' @return A string for the variance on the component
#' @seealso \code{\link[RGCCA]{rgcca}}, \code{\link[RGCCA]{sgcca}}
#' @examples
#' AVE = list(c(0.6, 0.5), c(0.7, 0.45))
#' rgcca_out = list(AVE = list(AVE_X = AVE))
#' # For the superblock (or the last block)
#' print_comp(rgcca_out, 1)
#' # "Axis 1 (70%)"
#' # For the first block
#' print_comp(rgcca_out, 2, 1)
#' # "Axis 2 (50%)"
#' @export
print_comp <- function(rgcca, n = NULL, i = NULL, outer = FALSE) {
    
    # by default, take the last block
    if (is.null(i))
        i <- length(rgcca$AVE$AVE_X)

    nvar <- sum(rgcca$a[[i]][, n] != 0)
    if (!is(rgcca, "sgcca") | nvar == length(rgcca$a[[i]][, n]))
        varText <- ""
    else
        varText <- paste0(nvar, " variables, ")
    
    ave <- quote(paste0(round(AVE[n] * 100, 1), "%"))
    if (isTRUE(outer)) {
        AVE <- rgcca$AVE$AVE_outer
        n <- c(1, 2)
        paste0("First outer comp. : ", paste(eval(ave), collapse = " & "))
    } else {
        AVE <- rgcca$AVE$AVE_X[[i]]
        paste0("Comp. ", n, " (", varText, eval(ave), ")")
    }
}
