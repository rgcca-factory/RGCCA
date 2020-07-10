#' Compute bootstrap (internal)
#'
#' Internal function for computing boostrap of RGCCA
#'
#' @inheritParams rgcca
#' @inheritParams plot_var_2D
#' @param rgcca_res Result of rgcca function
#' @return A list of RGCCA bootstrap weights
bootstrap_k <- function(
    rgcca_res,type="weight") {

    blocks_all <- rgcca_res$call$blocks
    rgcca_res_boot <- set_rgcca(rgcca_res, method="nipals",boot = TRUE)
    if(type=="weight")
    {
        w <- rgcca_res_boot$a        
    }
    if(type=="cor")
    {
       
                w=lapply(1:length(rgcca_res_boot$A),
                 function(j)
                {
                     
                     res=sapply(1:dim(rgcca_res_boot$A[[j]])[2],function(k){cor(rgcca_res_boot$Y[[j]][,1],rgcca_res_boot$A[[j]][,k],use="pairwise.complete.obs")})
                     names(res)=colnames(rgcca_res_boot$A[[j]])
                     return(res)
                 })
    }


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
    if(type=="weight")
    {
        check_sign_comp(rgcca_res, w)        
    }

}
