#' Correlation between the blocks
#' 
#' Calculates the pairwise correlation between each block for each of the 
#' components of the analysis
#'  
#' @inherit plot_ind
#' @inherit set_connection
#' @param comps A matrix containg the components of all the blocks
#' @examples
#' library(RGCCA)
#' data("Russett")
#' blocks = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
#'     politic = Russett[, 6:11] )
#' rgcca_out = rgcca.analyze(blocks)
#' get_cor_all(rgcca_out)
#' @return A list of matrix containg the pairwise correlation of the blocks
#' @export
get_cor_all <- function(
    rgcca, 
    blocks = rgcca$call$blocks, 
    comps = get_comp_all(rgcca)){

    comp <- list()

    for (i in seq(max(rgcca$call$ncomp))) {
        comp[[i]] <-  matrix(
            NA,
            NROW(comps),
            length(blocks),
            dimnames = list(rownames(comps), names(blocks))
        )

        for (n in names(blocks)) {
            pos <- grep(
                paste0(n, "_comp", i),
                colnames(comps))
            if (length(pos) > 0)
                comp[[i]][, n] <- comps[, pos]
        }

        comp[[i]] <- cor(comp[[i]], use = "pairwise.complete.obs")

    }

    return(comp)
}