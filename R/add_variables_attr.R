# Add removed variables in one submodel (bootstrap, cross-validation)
# to a list of weights and affect them to 0
add_variables_attr <- function(rgcca_res, w, type = "scale") {

    blocks_all <- list()
    for (x in seq(length(rgcca_res$call$blocks)))
        blocks_all[[x]] <- rgcca_res$call$blocks[[x]][intersect(rownames(rgcca_res$call$blocks[[1]]), rownames(w[[1]])), ]
    missing_var <- lapply(
        seq(length(w)), 
        function(x)
            setdiff(
                colnames(blocks_all[[x]]),
                names(w[[x]])))


    missing_tab <- lapply(
        seq(length(w)), 
        function(x) {
            #setNames(
                for (i in seq(length(missing_var[[x]])))
                    if (type == "scale")
                        0
                    else
                        mean(blocks_all[[x]][, i])
                #, missing_var[[x]])
        })
    
    for (i in seq(length(w))) {
        if (NROW(missing_tab[[i]]) != 0) {
            names(missing_tab[[i]]) <- missing_var[[i]]
             w[[i]] <- c(w[[i]], missing_tab[[i]])
        }
    }

    w <- lapply(seq(length(w)), function(x)
            w[[x]][colnames(blocks_all[[x]]), drop = FALSE])

    # names(w) <- names(rgcca_res$call$blocks)
    return(w)
}
