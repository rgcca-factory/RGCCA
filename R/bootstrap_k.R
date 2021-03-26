# Compute bootstrap (internal)
#
# Internal function for computing boostrap of RGCCA
#
# @param rgcca_res A fitted RGCCA object (see  \code{\link[RGCCA]{rgcca}})
# @return A list of RGCCA bootstrap weights/loadings.
bootstrap_k <- function(rgcca_res) {

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
                                          stop(error_message)
                                      }else{#Unknown message cases
                                      }
                                  })
        if (!(is.list(rgcca_res_boot))){
            nb_error = nb_error + 1
        }
    }
    
    if (is.list(rgcca_res_boot)){
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
    }else if (is.na(rgcca_res_boot)){
        return(NA)
    }
}