#' Get the components of the analysis
#' 
#' @inheritParams plot_ind
#' @inheritParams get_ctr
#' @param i_block_x An integer giving the index of a list of blocks
#' @param i_block_z An integer giving the index of a list of blocks (another
#' one, different from the one used in i_block_x)
#' @return A matrix containg each selected components and an associated response
get_comp <- function(
    rgcca_res,
    resp = rep(1, NROW(rgcca_res$Y[[1]])),
    compx = 1,
    compy = 2,
    compz = NULL,
    i_block_x = length(rgcca_res$Y),
    i_block_y = i_block_x,
    i_block_z = i_block_x,
    predicted = NULL){
    
    stopifnot(is(rgcca_res, "rgcca"))
    resp <- as.matrix(check_response(resp, rgcca_res$Y))

    for (i in c("i_block_x", "i_block_y", "i_block_z")) {
        if (!is.null(get(i)))
            check_blockx(i, get(i), rgcca_res$call$blocks)
    }
    check_ncol(rgcca_res$Y, i_block_x)
    for (i in c("x", "y", "z")) {
        j <- paste0("comp", i)
        if (!is.null(get(j)))
            check_compx(j, get(j), rgcca_res$call$ncomp, get(paste0("i_block_", i)))
    }

    df <- data.frame(
        rgcca_res$Y[[i_block_x]][, compx],
        rgcca_res$Y[[i_block_y]][, compy],
        rgcca_res$Y[[i_block_z]][, compz]
    )

    if (!is.null(predicted)) {
        df2 <- predicted[[2]][[i_block_x]][, c(compx, compy, compz)]
        colnames(df2) <- colnames(df)
        df1 <- df[rownames(df2),]
        df <- rbind(df1, df2)
        resp <- c(rep("obs",NROW(df2)),rep("pred",NROW(df2)))

    } else if (length(unique(resp)) > 1) {
        names <- row.names(resp)
        resp <- apply(as.matrix(resp), 1, as.character)

        if (!is.null(names)) {

            resp <- as.matrix(resp, row.names = names)
            name_blocks <- row.names(rgcca_res$Y[[i_block_x]])
            diff_column <- setdiff(name_blocks, names)

            if (identical(diff_column, name_blocks)) {
                warning("No match has been found with the row names of the group file.")
                resp <- rep("NA", NROW(df))

            } else {
                if (length(diff_column) > 0) {
                    resp[diff_column] <- "NA"
                    names(resp)[names(resp) == ""] <- names
                } else {
                    names(resp) <- names
                }
                resp <- resp[row.names(rgcca_res$Y[[i_block_x]])]
            }

        } else {
            # warning("No rownames have been found in the group file. The rownames of the selected block of RGCCA have been used.")
            # resp <- rep("NA", NROW(df))
            # rownames(resp) <- rownames(rgcca_res$A[[i_block_x]])
            if (length(resp) != NROW(rgcca_res$A[[i_block_x]]))
                stop("resp argument should have the same size than the number of rows in the selected block.")
        }
    } else
        resp <- resp[seq(NROW(df)), ]

    if ((!is.character2(resp) &&
        length(unique(resp)) > 5) || 
            unique(resp) == 1 ) {
        resp[resp == "NA"] <- NA
        df$resp <- as.numeric(resp)
    }else
        df$resp <- as.factor(resp)

    return(df)
}
