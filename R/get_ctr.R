#' Variable contribution
#' 
#' Extract the contibution of variables to the model by using correlation or weight
#' 
#' @inheritParams plot_var_2D
#' @param type A character giving the choice ot the index between cor or weight
#' @param compz An integer giving the index of the analysis component used
#' for the z-axis
#' @param collapse_each A bolean, if collapse is actived, to get the correlation 
#' of each blocks to the i_block (FALSE, default) or to itself (TRUE)
#' @return A dataframe containing the indexes for each selected components
get_ctr <- function(
    rgcca_res,
    compx = 1,
    compy = 2,
    compz = NULL,
    i_block = length(rgcca_res$call$blocks),
    type = "cor",
    collapse = FALSE,
    collapse_auto = FALSE) {

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
        if (collapse_auto)
            i_block <- y
        cor(
            blocks[[y]][rownames(rgcca_res$Y[[y]]), ],
            rgcca_res$Y[[i_block]][, x],
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
