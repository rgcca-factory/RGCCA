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
    rgcca,
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
    response = NULL) {

    if (is.null(blocks)) {
        blocks <- intersection(rgcca$call$blocks) -> blocks.all
        connection <- rgcca$call$connection
        ncomp <- rgcca$call$ncomp
        scheme <- rgcca$call$scheme
        bias <- rgcca$call$bias
        tol <- rgcca$call$tol
        superblock <- rgcca$call$superblock
        type <- rgcca$call$type
        init <- rgcca$call$init
        sameBlockWeight <- rgcca$call$scale
        scale <- rgcca$call$scale

        if (rgcca$call$type %in% c("sgcca","spls","spca")) {
            penalty <- rgcca$call$c1
            par <- "sparsity"
        } else {
            penalty <- rgcca$call$tau
            par <- "tau"
        }

        if (superblock) {
            blocks <- blocks[-length(blocks)]
            connection <- NULL
        }

        if (!is.null(rgcca$call$response))
            response <- length(rgcca$call$blocks)

    } else {
        blocks.all <- intersection(blocks)

        if (tolower(type) %in% c("sgcca", "spca", "spls")) {
            if (!missing(tau))
               stop(paste0("penalty parameter required for ", tolower(type), "."))
            par <- "sparsity"
            penalty <- sparsity
        } else {
            if (!missing(sparsity))
               stop(paste0("tau parameter required for ", tolower(type), "."))
            par <- "tau"
            penalty <- tau
        }
    }

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

    w <- check_sign_comp(rgcca, w)

    return(w)
}
