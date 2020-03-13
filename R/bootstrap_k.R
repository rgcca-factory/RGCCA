#' Compute bootstrap (internal)
#'
#' Internal function for computing boostrap of RGCCA
#'
#' @inheritParams rgcca
#' @inheritParams plot_var_2D
#' @param rgcca_res Result of rgcca function
#' @return A list of RGCCA bootstrap weights
bootstrap_k <- function(
    rgcca_res,
    ...) {

    blocks_all <- rgcca_res$call$blocks
    rgcca_res <- set_rgcca(rgcca_res, boot = TRUE, ...)
    w <- rgcca_res$a

    # Add removed variables
    missing_var <- lapply(
        seq(length(w)),
        function(x)
            setdiff(colnames(blocks_all[[x]]), rownames(w[[x]]))
    )

    missing_tab <- lapply(
        seq(length(w)),
        function(x)
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

    w <- lapply(
        seq(length(w)), 
        function(x) 
            w[[x]][colnames(blocks_all[[x]]),
            ,
            drop = FALSE])

    names(w) <- names(rgcca_res$call$blocks)

    check_sign_comp(rgcca_res, w)
}
