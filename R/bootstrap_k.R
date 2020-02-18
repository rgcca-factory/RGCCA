#' Compute bootstrap (internal)
#'
#' Internal function for computing boostrap of RGCCA
#'
#' @inheritParams rgcca
#' @inheritParams plot_var_2D
#' @return A list of RGCCA bootstrap weights
#' @examples
#' library(RGCCA)
#' data("Russett")
#' blocks = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
#'     politic = Russett[, 6:11] )
#' rgcca_out = rgcca(blocks)
#' bootstrap_k(rgcca_out)
#' bootstrap_k(rgcca_out, lapply(blocks, scale), superblock = FALSE)
#' @export
bootstrap_k <- function(
    rgcca_res,
    blocks = NULL,
    connection = 1 - diag(length(blocks)),
    tau = rep(1, length(blocks)),
    sparsity = rep(1, length(blocks)),
    ncomp = rep(2, length(blocks)),
    scheme = "factorial",
    init = "svd",
    bias = TRUE,
    tol = 1e-03,
    type = "rgcca",
    scale = TRUE,
    sameBlockWeight = TRUE,
    superblock = TRUE,
    response = NULL,
    method = "complete") {

    if (is.null(blocks)) {
        blocks <- rgcca_res$call$blocks -> blocks.all
        method <- rgcca_res$call$method
        connection <- rgcca_res$call$connection
        ncomp <- rgcca_res$call$ncomp
        scheme <- rgcca_res$call$scheme
        bias <- rgcca_res$call$bias
        tol <- rgcca_res$call$tol
        superblock <- rgcca_res$call$superblock
        type <- rgcca_res$call$type
        init <- rgcca_res$call$init
        sameBlockWeight <- rgcca_res$call$sameBlockWeight 
        scale <- rgcca_res$call$scale

        if (rgcca_res$call$type %in% c("sgcca","spls","spca")) {
            penalty <- rgcca_res$call$sparsity
            par <- "sparsity"
        } else {
            penalty <- rgcca_res$call$tau
            par <- "tau"
        }

        if (superblock) {
            blocks <- blocks[-length(blocks)]
            connection <- NULL
        }

        if (!is.null(rgcca_res$call$response))
            response <- length(rgcca_res$call$blocks)

    } else {
        blocks.all <- blocks

        if (tolower(type) %in% c("sgcca", "spca", "spls")) {
            if (!missing(tau) & missing(sparsity)){stop(paste0("sparsity parameter required for ", tolower(type), "instead of tau."))}
               
            par <- "sparsity"
            penalty <- sparsity
        } else {
            if (!missing(sparsity) & missing(tau))
            {
                stop(paste0("tau parameter required for ", tolower(type), "instead of sparsity."))
            }
              
            par <- "tau"
            penalty <- tau
        }
    }
print(penalty)
    boot_blocks <- list(NULL)
    while (any(sapply(boot_blocks, function(x) length(x)) == 0)) {

        id_boot <- sample(NROW(blocks[[1]]), replace = TRUE)

        boot_blocks <- lapply(
            blocks, 
            function(x) x[id_boot, , drop = FALSE])

        boot_blocks <- remove_null_sd(boot_blocks)
    }

    func <- quote(
        rgcca(
            boot_blocks,
            connection,
            superblock = superblock,
            response = response,
            ncomp = ncomp,
            scheme = scheme,
            scale = scale,
            sameBlockWeight = sameBlockWeight,
            type = type,
            verbose = FALSE,
            init = init,
            bias = bias,
            method=method,
        tol = tol
    ))
    
    func[[par]] <- penalty
    w <- eval(as.call(func))$a

    # Add removed variables
    missing_var <- lapply(
            seq(length(w)),
            function(x) setdiff(colnames(blocks.all[[x]]), rownames(w[[x]]))
        )

    missing_tab <- lapply(
        seq(length(w)),
        function(x)
            matrix(
                0,
                length(missing_var[[x]]),
                ncomp[x],
                dimnames = list(missing_var[[x]], seq(ncomp[x]))
        ))

    # bug mapply with pca
    w <- lapply(seq(length(missing_tab)), function(x) {
        if (NROW(missing_tab[[x]]) != 0)
            rbind(w[[x]], missing_tab[[x]])
        else
            w[[x]]
        })

    w <- lapply(seq(length(w)), function(x) w[[x]][colnames(blocks.all[[x]]), , drop = FALSE])

    names(w) <- names(blocks.all)

    w <- check_sign_comp(rgcca_res, w)

    return(w)
}
