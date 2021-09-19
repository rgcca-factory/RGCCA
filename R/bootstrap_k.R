# Compute bootstrap (internal)
#
# Internal function for computing boostrap of RGCCA
#
# @param rgcca_res A fitted RGCCA object (see  \code{\link[RGCCA]{rgcca}})
# @return A list of RGCCA bootstrap weights/loadings.
bootstrap_k <- function(rgcca_res, boot_blocks = NULL, column_sd_null = NULL) {
    rgcca_res_boot <- set_rgcca(rgcca_res,
                                blocks    = boot_blocks,
                                NA_method = "nipals",
                                boot      =  is.null(boot_blocks))
    #block-weight vector
    W = add_variables_submodel(rgcca_res, rgcca_res_boot$a)

    #block-loadings vector
    A = check_sign_comp(rgcca_res, rgcca_res_boot$a)

    Y = lapply(seq_along(W),
               function(j) pm(rgcca_res_boot$call$blocks[[j]], A[[j]])
               )
    L <- lapply(seq_along(W),
                function(j)
                    cor(rgcca_res_boot$call$blocks[[j]], Y[[j]],
                        use = "pairwise.complete.obs")
                )

    names(L) = names(rgcca_res$a)
    return(list(W = W, L = L))
}
