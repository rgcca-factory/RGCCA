# Compute bootstrap (internal)
#
# Internal function for computing boostrap of RGCCA
#
# @inheritParams rgcca
# @inheritParams plot_var_2D
# @return A list of RGCCA bootstrap weights
bootstrap_k <- function(rgcca_res, type = "weight") {

    rgcca_res_boot = NA
    nb_error       = 0
    while( !(is.list(rgcca_res_boot)) && (nb_error < 5)){
        rgcca_res_boot = tryCatch(set_rgcca(rgcca_res, NA_method = "nipals", boot = TRUE),
                                  error = function(error_message){
                                      if (grepl("L1/L2 projection issue", toString(error_message), fixed = TRUE)){
                                          if(nb_error == 4){
                                              warning(error_message)
                                          }
                                          return(NA)
                                      }else{#Unknown message cases
                                          stop(error_message)
                                      }
                                  })
        if (!(is.list(rgcca_res_boot))){
            nb_error = nb_error + 1
        }
    }
    
    if (is.list(rgcca_res_boot)){
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
    }else if (is.na(rgcca_res_boot)){
        return(NA)
    }
}