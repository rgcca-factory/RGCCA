#' Variable contribution
#' 
#' Extract the contibution of variables to the model by using correlation or weight
#' 
#' @inheritParams plot_var_2D
#' @param type A character giving the choice ot the index between cor or weight
#' @param compz An integer giving the index of the analysis component used
#' for the z-axis
#' @examples
#' library(RGCCA)
#' data("Russett")
#' blocks = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
#'     politic = Russett[, 6:11] )
#' rgcca_out = rgcca.analyze(blocks, ncomp = c(3,2,4))
#' get_ctr(rgcca_out)
#' # On the first block and with weights
#' get_ctr(rgcca_out, 2, 1, i_block = 1, type = "weight")
#' # With 3 components and on the variables of two blocks
#' superblocks <- rep(list(Reduce(cbind, c(blocks[1], blocks[3]))), 2)
#' names(superblocks) <- names(blocks)[c(1, 3)]
#' rgcca_out = rgcca.analyze(blocks[c(1,3)], ncomp = c(3,4))
#' rgcca_out$blocks = superblocks
#' get_ctr(rgcca_out, compz = 3, i_block = 1, type = "cor", collapse = TRUE)
#' get_ctr(rgcca_out, 2, 1, 3, 1, "weight", TRUE)
#' @return A dataframe containing the indexes for each selected components
#' @export
get_ctr <- function(
    rgcca,
    compx = 1,
    compy = 2,
    compz = NULL,
    i_block = length(rgcca$blocks),
    type = "cor",
    collapse = FALSE) {
    
    match.arg(type, c("cor", "weight"))
    stopifnot(!missing(rgcca))

    if (!collapse)
        i_block_2 <- i_block
    else
        i_block_2 <- 1

    row.names = colnames(rgcca$blocks[[i_block]])

    if (type == "cor")
        f <- function(x){
            cor(
                rgcca$blocks[[i_block_2]],
                rgcca$Y[[i_block]][, x],
                use = "pairwise.complete.obs"
            )
        }
    else{
        if (!collapse)
            f <- function(x) rgcca$a[[i_block]][, x]
        else
            f <- function(x) unlist(
                sapply(
                    seq(length(rgcca$blocks)),
                    function(y) rgcca$a[[y]][, x]
                )
            )
    }

    data.frame(
        sapply(
            c(compx, compy, compz[compz >= rgcca$ncomp[i_block]]),
            function(x) f(x)
        ),
        row.names = row.names
    )

}
