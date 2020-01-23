#' Variable contribution
#' 
#' Extract the contibution of variables to the model by using correlation or weight
#' 
#' @inheritParams plot_var_2D
#' @param type A character giving the choice ot the index between cor or weight
#' @param compz An integer giving the index of the analysis component used
#' for the z-axis
#' @param i_block_2 An integer giving the index of a list of blocks to be 
#' correlated to i_block if this option is selected (default to i_block)
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
#' rgcca_out$call$blocks = superblocks
#' get_ctr(rgcca_out, compz = 3, i_block = 1, type = "cor", collapse = TRUE)
#' get_ctr(rgcca_out, 2, 1, 3, 1, "weight", TRUE)
#' @return A dataframe containing the indexes for each selected components
#' @export
get_ctr <- function(
    rgcca,
    compx = 1,
    compy = 2,
    compz = NULL,
    i_block = length(rgcca$call$blocks),
    type = "cor",
    collapse = FALSE,
    i_block_2 = NULL) {

    match.arg(type, c("cor", "weight"))
    stopifnot(!missing(rgcca))

    blocks <- rgcca$call$blocks

    if (!collapse) {
        row.names <- colnames(blocks[[i_block]])
    }else{
        if (rgcca$call$superblock)
            blocks <- blocks[-length(blocks)]
        row.names <- unlist(lapply(blocks, colnames))
    }

    if (type == "cor") {
        if (!collapse)
            f <- function(x){
                    if (is.null(i_block_2))
                        i_block_2 <- i_block
                    cor(
                        blocks[[i_block_2]][rownames(rgcca$Y[[i_block]]), ],
                        rgcca$Y[[i_block]][, x],
                        use = "pairwise.complete.obs"
                    )
            }
        else
            f <- function(x){
                unlist(
                    lapply(
                        seq(length(blocks)),
                        function(y){    
                            if (is.null(i_block_2))
                                i_block_2 <- y
                            cor(
                                blocks[[i_block_2]][rownames(rgcca$Y[[y]]), ],
                                rgcca$Y[[y]][, x],
                                use = "pairwise.complete.obs"
                            )
                        }
                    )
                )
            }
    }else{
        if (!collapse)
            f <- function(x) rgcca$a[[i_block]][, x]
        else
            f <- function(x) unlist(
                lapply(
                    seq(length(blocks)),
                    function(y) rgcca$a[[y]][, x]
                )
            )
    }

    res <- data.frame(
        sapply(
            c(compx, compy, compz),
            function(x) f(x),
            simplify = FALSE
        ),
        row.names = row.names
    )
    colnames(res) <- seq(NCOL(res)) # for save_var
    return(res)

}

#[compz >= rgcca$call$ncomp[i_block]]