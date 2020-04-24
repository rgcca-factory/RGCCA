set_rgcca <- function(
    rgcca_res,
    blocks = NULL,
    connection = NULL,
    tau = 1,
    sparsity = 1,
    ncomp = 2,
    scheme = "factorial",
    init = "svd",
    bias = TRUE,
    tol = 1e-03,
    type = "rgcca",
    scale = TRUE,
    sameBlockWeight = TRUE,
    superblock = FALSE,
    response = NULL,
    method = "complete",
    boot = FALSE,
    inds = NULL) {

    if (is.null(blocks)) {
        blocks <- rgcca_res$call$blocks
        method <- rgcca_res$call$method
        connection <- rgcca_res$call$connection
        scheme <- rgcca_res$call$scheme
        bias <- rgcca_res$call$bias
        superblock <- rgcca_res$call$superblock
        type <- rgcca_res$call$type
        init <- rgcca_res$call$init
        scale <- FALSE
        sameBlockWeight <- FALSE

        if (superblock) {
            J <- length(blocks)
            blocks <- blocks[-J]
            for (i in c("tau", "sparsity", "ncomp")) {
                if (class(rgcca_res$call[[i]]) %in% c("matrix", "data.frame"))
                  rgcca_res$call[[i]] <- rgcca_res$call[[i]][,-J]
                else
                  rgcca_res$call[[i]] <- rgcca_res$call[[i]][-J]
            }

            connection <- NULL
        }

        sparsity <- rgcca_res$call$sparsity
        tau <- rgcca_res$call$tau
        ncomp <- rgcca_res$call$ncomp

        if (!is.null(rgcca_res$call$response))
            response <- length(rgcca_res$call$blocks)

    }else
        blocks <- scaling(blocks, scale, sameBlockWeight = sameBlockWeight)

    if (!boot)
        blocks <- intersection(blocks)

    if (tolower(type) %in% c("sgcca", "spca", "spls")) {

        if (!is.null(blocks) && !missing(tau) && missing(sparsity))
            stop(paste0("sparsity parameter required for ", tolower(type), "instead of tau."))

        par <- "sparsity"
        penalty <- sparsity

    } else {

        if (!is.null(blocks) && !missing(sparsity) && missing(tau))
            stop(paste0("tau parameter required for ", tolower(type), "instead of sparsity."))

        par <- "tau"
        penalty <- tau
    }

    if (boot) {
        boot_blocks <- list(NULL)
        while (any(sapply(boot_blocks, function(x) length(x)) == 0)) {

            id_boot <- sample(NROW(blocks[[1]]), replace = TRUE)

            boot_blocks <- lapply(
                blocks, 
                function(x) x[id_boot, , drop = FALSE])

            boot_blocks <- remove_null_sd(boot_blocks)
        }
    }else
        boot_blocks <- lapply(blocks, function(x) x[-inds, , drop = FALSE])

    func <- quote(
        rgcca(
            boot_blocks,
            connection,
            superblock = superblock,
            response = response,
            ncomp = ncomp,
            scheme = scheme,
            scale = FALSE,
            sameBlockWeight = FALSE,
            type = type,
            verbose = FALSE,
            init = init,
            bias = bias,
            method = method,
            tol = tol
        ))

    func[[par]] <- penalty

    res <- eval(as.call(func))
    attributes(res)$bigA_scaled <- blocks
    return(res)
}
