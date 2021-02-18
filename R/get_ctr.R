# Variable contribution
# 
# Extract the contibution of variables to the model by using correlation or weight
# 
# @inheritParams plot_var_2D
# @inheritParams plot_var_1D
# @param compz An integer giving the index of the analysis component used
# for the z-axis
# @param i_block_2 An integer giving the index of a list of blocks to be 
# correlated to i_block if this option is selected
# @return A dataframe containing the indexes for each selected components
get_ctr <- function(
    rgcca_res,
    compx = 1,
    compy = 2,
    compz = NULL,
    i_block = length(rgcca_res$blocks),
    type = "cor",
    collapse = FALSE,
    i_block_2 = i_block) {

    match.arg(type, c("cor", "weight"))
    stopifnot(!missing(rgcca_res))

    blocks <- rgcca_res$blocks
    y <- NULL

    if (collapse && rgcca_res$call$superblock) {
        blocks <- blocks[-length(blocks)]
    }

    if (type == "cor")
        f2 <- function(x, y) cor2(blocks[[y]], rgcca_res$Y[[i_block]][, x], 
                                  use = "pairwise.complete.obs")
    else
        f2 <- function(x, y) rgcca_res$a[[y]][, x]

    if (!collapse)
        f <- function(x)
            f2(x, i_block_2)
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
                    stop_rgcca("The index of the selected analysis component doesn't exist.")
                f(x)
            },
            simplify = FALSE
        )
    )
    colnames(res) <- seq(NCOL(res)) # for save_var
    return(res)

}
