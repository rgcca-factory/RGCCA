# Compute bootstrap (internal)
#
# Internal function for computing boostrap of RGCCA
#
# @param rgcca_res A fitted RGCCA object (see  \code{\link[RGCCA]{rgcca}})
# @return A list of RGCCA bootstrap weights/loadings.
bootstrap_k <- function(rgcca_res) {
    rgcca_res_boot <- set_rgcca(rgcca_res, NA_method = "nipals", boot = TRUE)

    #block-weight vector
    W = add_variables_submodel(rgcca_res, rgcca_res_boot$a)
    boot_data = add_variables_data(rgcca_res, rgcca_res_boot$call$blocks, boot = TRUE)

    Y = lapply(seq_along(W), function(j) pm(boot_data[[j]], W[[j]]))

    options(warn = -1)
    L <- lapply(
        seq_along(W),
        function(j) {
            y <- cor(boot_data[[j]], Y[[j]],
                use = "pairwise.complete.obs")
            y[is.na(y)] <- 0
            y
        }
    )
    options(warn = 0)

    names(L) = names(rgcca_res$a)
    return(list(W = W, L = L))
}
