# An intern function used by sgcca.permute to perform multiple sgcca with permuted rows
# data("Russett")
# blocks = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
#     politic = Russett[, 6:11] )
# rgcca_out = rgcca(blocks)
# rgcca_permutation_k(blocks)
rgcca_permutation_k <- function(
    blocks,
    par = list("ncomps", expand.grid(rep(list(seq(2)), length(blocks)))),
    connection = 1 - diag(length(blocks)),
    response = NULL,
    tau = rep(1, length(blocks)),
    sparsity = rep(1, length(blocks)),
    ncomp = rep(2, length(blocks)),
    scheme = "factorial",
    scale = TRUE,
    init = "svd",
    bias = TRUE,
    tol = 1e-03,
    type = "rgcca",
    superblock = TRUE,
    perm = TRUE,
    n_cores = parallel::detectCores() - 1) {

    if (perm) {
        for (k in seq(length(blocks)))
            blocks[[k]] <- as.matrix(blocks[[k]][sample(seq(NROW(blocks[[k]]))), ])
    }

    simplify2array(
        parallel::mclapply(
            seq(NROW(par[[2]])),
            function(i) {
                func <- quote(
                    rgcca(
                        blocks = blocks,
                        connection = connection,
                        response = response,
                        superblock = superblock,
                        ncomp = ncomp,
                        scheme = scheme,
                        scale = scale,
                        type = type,
                        init = init,
                        bias = bias,
                        tol = tol,
                        quiet = TRUE,
                        method = "complete"
                    ))
                switch(
                    par[[1]],
                    "ncomp" = {
                        ncomp <- par[[2]][i, ]
                            if (tolower(type) %in% c("sgcca", "spca", "spls"))
                                func[["penalty"]] <- sparsity
                            else
                                func[["tau"]] <- tau
                    },
                    "sparsity" = {
                        func[[par[[1]]]] <- par[[2]][i, ]
                        type <- "sgcca"
                })
                crit <- eval(as.call(func))$crit
                
                return(sum(sapply(crit, function(x) sum(x))))
            },
        mc.cores = n_cores))
}
