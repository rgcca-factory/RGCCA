#' Get the indexes of the analysis
#' 
#' @inheritParams plot_var_2D
#' @inheritParams get_ctr
#' @return A matrix containg the indexes (correlation of the blocks with a 
#' component or their weights) for each selected component and an associated response

get_ctr2 <- function(
    rgcca_res,
    compx = 1,
    compy = 2,
    compz = NULL,
    i_block = length(rgcca_res$call$blocks),
    type = "cor",
    n_mark = 100,
    collapse = FALSE,
    remove_var = TRUE,
    resp=NULL) {

    stopifnot(is(rgcca_res, "rgcca"))
    check_blockx("i_block", i_block, rgcca_res$call$blocks)
    check_ncol(rgcca_res$a, i_block)
    for (i in c("compx", "compy", "compz")) {
        if (!is.null(get(i)))
            check_compx(i, get(i), rgcca_res$call$ncomp, i_block)
    }
    check_integer("n_mark", n_mark)
    for (i in c("collapse", "remove_var"))
        check_boolean(i, get(i))

    x <- y <- selectedVar <- NULL
    blocks <- rgcca_res$call$blocks

    if (collapse) {
        if (rgcca_res$call$superblock) {

            blocks <- blocks[-length(blocks), drop = FALSE]
            if (i_block > length(blocks))
                i_block <- length(blocks)
        }
        rgcca_res$call$superblock <- TRUE
        blocks.all <- blocks
        blocks <- rep(list(Reduce(cbind, blocks)), length(blocks))
        names(blocks) <- names(blocks.all)
    }

    df <- get_ctr(rgcca_res, compx, compy, compz, i_block, type, collapse)

    if (rgcca_res$call$type %in% c("spls", "spca", "sgcca") ){

        if (collapse)
            J <- seq(length(rgcca_res$a))
        else
            J <- i_block

        selectedVar <- unlist(
            lapply(
                J,
                function(x)
                    apply(
                        sapply(
                            c(compx, compy, compz[compz >= rgcca_res$call$ncomp[x]]),
                            function(y) rgcca_res$a[[x]][, y] != 0),
                        1,
                        function(z) Reduce("|", z)
                    )
                )
            )
        df <- df[ names(which(selectedVar)), ]

    }

    if (n_mark > NROW(df))
        n_mark <- NROW(df)

    # TODO: function in other place
    if (remove_var) {
        selectedVar <- unique(as.vector(
            sapply(seq(length(c(compx, compy, compz))), function(x)
                row.names(data.frame(df[order(abs(df[, x]), decreasing = TRUE), ])[seq(n_mark), ]))
        ))
        df <- df[selectedVar, ]
    } else
        selectedVar <- row.names(df)

    # group by blocks
    if(is.null(resp))
    {
        if (rgcca_res$call$superblock & (collapse | (i_block == length(rgcca_res$a)))) 
        {
            if (collapse)
                resp <- get_bloc_var(lapply(blocks.all, t), TRUE)
            else{
                resp <- get_bloc_var(rgcca_res$a)
                
                resp <- resp[
                    unlist(
                        lapply(
                            seq(length(selectedVar)),
                            function(x) which(colnames(blocks[[length(blocks)]]) == selectedVar[x])
                        )
                    )
                    ]
            }
            # df <- resp[row.names(df)]
            
        } else
            resp <- rep(1, NROW(df))
        
    }
   
    data.frame(df, resp)
}
