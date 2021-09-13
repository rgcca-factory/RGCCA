# Add removed variables in one submodel (bootstrap, cross-validation)
# to a list of weights and affect them to 0
add_variables_submodel <- function(rgcca_res, w) {

    blocks_all <- rgcca_res$call$blocks
    missing_var <- lapply(seq(length(w)), function(x)
            setdiff(colnames(blocks_all[[x]]),
                rownames(w[[x]])))

    missing_tab <- lapply(seq(length(w)), function(x)
            matrix(
                0,
                length(missing_var[[x]]),
                rgcca_res$call$ncomp[x],
                dimnames = list(missing_var[[x]], seq(rgcca_res$call$ncomp[x]))
            ))

    # bug mapply with pca
    w <- lapply(seq(length(missing_tab)), function(x) {
        if (NROW(missing_tab[[x]]) != 0)
            rbind(w[[x]], missing_tab[[x]])
        else
            w[[x]]
    })

    w <- lapply(seq(length(w)), function(x)
            w[[x]][colnames(blocks_all[[x]]), , drop = FALSE])

    names(w) <- names(rgcca_res$call$blocks)
    return(check_sign_comp(rgcca_res, w))
}
