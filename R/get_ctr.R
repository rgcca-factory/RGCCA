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
#' rgcca_out = rgcca(blocks, ncomp = c(3,2,4))
#' get_ctr(rgcca_out)
#' # On the first block and with weights
#' get_ctr(rgcca_out, 2, 1, i_block = 1, type = "weight")
#' # With 3 components and on the variables of two blocks
#' rgcca_out = rgcca(blocks[c(1,3)], ncomp = c(3,4))
#' get_ctr(rgcca_out, compz = 3, i_block = 1, type = "cor", collapse = TRUE)
#' get_ctr(rgcca_out, 2, 1, 3, 1, "weight", TRUE)
#' @return A dataframe containing the indexes for each selected components
#' @export
get_ctr <- function(
    rgcca_res,
    compx = 1,
    compy = 2,
    compz = NULL,
    i_block = length(rgcca_res$call$blocks),
    type = "cor",
    collapse = FALSE,
    i_block_2 = NULL) {

    match.arg(type, c("cor", "weight"))
    stopifnot(!missing(rgcca_res))

    blocks <- rgcca_res$call$blocks
    y <- NULL

    if (!collapse)
        row.names <- colnames(blocks[[i_block]])
    else{
        if (rgcca_res$call$superblock)
            blocks <- blocks[-length(blocks)]
        row.names <- unlist(lapply(blocks, colnames))
    }

    if (type == "cor")
        f2 <- function(x, y){    
        if (is.null(i_block_2))
            i_block_2 <- y
        cor(
            blocks[[i_block_2]][rownames(rgcca_res$Y[[y]]), ],
            rgcca_res$Y[[y]][, x],
            use = "pairwise.complete.obs"
        )
    }
    else
        f2 <- function(x, y) rgcca_res$a[[y]][, x]

    if (!collapse)
        f <- function(x)
            f2(x, i_block)
    else
        f <- function(x){
            unlist(
                lapply(
                    seq(length(blocks)),
                    function(y) f2(x, y)
                )
            )
        }

    res <- data.frame(
        sapply(
            c(compx, compy, compz),
            function(x){
                if (x > rgcca_res$call$ncomp[i_block])
                    stop("The index of the selected analysis component doesn't exist.")
                f(x)
            },
            simplify = FALSE
        ),
        row.names = row.names
    )
    colnames(res) <- seq(NCOL(res)) # for save_var
    return(res)

}
