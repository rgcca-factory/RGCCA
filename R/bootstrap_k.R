# Compute bootstrap (internal)
#
# Internal function for computing boostrap of RGCCA
#
# @inheritParams rgcca
# @inheritParams plot_var_2D
# @return A list of RGCCA bootstrap weights
bootstrap_k <- function(rgcca_res, type = "weight") {

    rgcca_res_boot <- set_rgcca(rgcca_res, method = "nipals", boot = TRUE)

    if (type == "weight") {
        add_variables_submodel(rgcca_res, rgcca_res_boot$a)

    } else {
        w <- lapply(1:length(rgcca_res_boot$call$blocks), function(j) {
            res <- sapply(1:dim(rgcca_res_boot$call$blocks[[j]])[2], function(k) {
                cor(
                    rgcca_res_boot$Y[[j]][, 1],
                    rgcca_res_boot$call$blocks[[j]][, k],
                    use = "pairwise.complete.obs")
            })
            names(res) <- colnames(rgcca_res_boot$call$blocks[[j]])
            return(res)
        })
    }
}
