#' Compute bootstrap (internal)
#'
#' Internal function for computing boostrap of RGCCA
#'
#' @inheritParams rgcca.analyze
#' @inheritParams plot_var_2D
#' @return A list of RGCCA bootstrap weights
#' @examples
#' library(RGCCA)
#' data("Russett")
#' blocks = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
#'     politic = Russett[, 6:11] )
#' rgcca_out = rgcca.analyze(blocks)
#' bootstrap_k(rgcca_out)
#' bootstrap_k(rgcca_out, lapply(blocks, scale), superblock = FALSE)
#' @export
bootstrap_k <- function(
    rgcca,
    blocks = NULL,
    connection = 1 - diag(length(blocks)),
    tau = rep(1, length(blocks)),
    ncomp = rep(2, length(blocks)),
    scheme = "factorial",
    init = "svd",
    bias = TRUE,
    tol = 1e-03,
    type = "rgcca",
    superblock = TRUE) {

    if (is.null(blocks))
        blocks.all <- rgcca$call$blocks
    else
        blocks.all <- blocks

    if (is.null(blocks)) {
        blocks <- rgcca$call$blocks
        connection <- rgcca$call$connection
        ncomp <- rgcca$call$ncomp
        scheme <- rgcca$call$scheme
        bias <- rgcca$call$bias
        tol <- rgcca$call$tol
        superblock <- rgcca$call$superblock
        type <- class(rgcca)
        init <- rgcca$init

        if (class(rgcca) %in% c("sgcca", "spls", "spca"))
            tau <- rgcca$call$c1
        else
            tau <- rgcca$call$tau

        if (superblock) {
            blocks <- blocks[-length(blocks)]
            connection <- NULL
        }
    }

    # Shuffle rows
    id_boot <- sample(NROW(blocks[[1]]), replace = TRUE)

    if (any(sapply(blocks, function(x) is.null(attr(x, 'scaled:center')))))
            stop("Blocks should be scaled before performing bootstraps.")
    else
        boot_blocks <- lapply(
            blocks, 
            function(x) scale2(x[id_boot, , drop = FALSE], scale = FALSE))

    boot_blocks <- remove_null_sd(boot_blocks)

    # Get boostraped weights
    w <- rgcca.analyze(
        boot_blocks,
        connection,
        superblock = superblock,
        tau = tau,
        ncomp = ncomp,
        scheme = scheme,
        scale = FALSE,
        type = type,
        verbose = FALSE,
        init = init,
        bias = bias,
        tol = tol
    )$a

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
