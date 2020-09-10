# Add removed variables in one submodel (bootstrap, cross-validation)
# to a list of weights and affect them to 0
add_variables_data <- function(rgcca_res, w) {

    blocks_all <- list()
    for (x in seq(length(rgcca_res$call$blocks)))
    {
        blocks_all[[x]] <- rgcca_res$call$blocks[[x]][intersect(rownames(rgcca_res$call$blocks[[1]]), rownames(w[[1]])), ]
        if(is.null(dim(blocks_all[[x]])))
        {
            b=blocks_all[[x]]
            blocks_all[[x]]=matrix(blocks_all[[x]],ncol=1)
            rownames(blocks_all[[x]])=names(b)
            colnames(blocks_all[[x]])=names(rgcca_res$call$blocks)[x]
        }
    }
     
    missing_var <- lapply(
        seq(length(w)), function(x)
            setdiff(
                colnames(blocks_all[[x]]),
                colnames(w[[x]])))

    missing_tab <- lapply(seq(length(w)), function(x)
    {
        if(!is.null(dim(x))&&dim(x)[2]==0)
        {
            M=matrix(
                0,
                nrow(blocks_all[[x]]),
                length(missing_var[[x]]),
                dimnames = list(rownames(blocks_all[[x]]), missing_var[[x]])
            ) 
        }else
        {
            M=matrix(
                0,
                nrow(blocks_all[[x]]),
                length(missing_var[[x]]),
                dimnames = list(names(blocks_all[[x]]), missing_var[[x]])
            )
        }
      
        return(M)  
    }
    )

    # bug mapply with pca

    w <- lapply(seq(length(missing_tab)), function(x) {
        if (NCOL(missing_tab[[x]]) != 0)
            cbind(w[[x]], missing_tab[[x]])
        else
            w[[x]]
    })

    w <- lapply(seq(length(w)), function(x)
            w[[x]][, colnames(blocks_all[[x]]), drop = FALSE])
    names(w) <- names(rgcca_res$call$blocks)
    return(w)
}
