# Compute bootstrap (internal)
#
# Internal function for computing boostrap of RGCCA
#
# @param rgcca_res A fitted RGCCA object (see  \code{\link[RGCCA]{rgcca}})
# @return A list of RGCCA bootstrap weights/loadings.
bootstrap_k <- function(rgcca_res) {
    rgcca_res_boot <- set_rgcca(rgcca_res, NA_method = "nipals", boot = TRUE)
    W = add_variables_submodel(rgcca_res, rgcca_res_boot$a)
    L <- lapply(seq_along(rgcca_res$a), 
                function(x) 
                    cor(rgcca_res_boot$call$blocks[[x]], 
                        rgcca_res_boot$Y[[x]], 
                        use = "pairwise.complete.obs")
    )
    names(L) = names(rgcca_res$a)
    L = check_sign_comp(rgcca_res, L)
    return(list(W = W, L = L))
}