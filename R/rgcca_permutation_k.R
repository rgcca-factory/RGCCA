# An intern function used by sgcca.permute to perform multiple sgcca with permuted rows
# data("Russett")
# blocks = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
#     politic = Russett[, 6:11] )
# rgcca_permutation_k(blocks)
rgcca_permutation_k <- function(
    blocks,
    par = "ncomp",
    par_value=rep(1,length(blocks)),
    tol = 1e-03,
    type = "rgcca",
    sparsity = rep(1, length(blocks)),
    perm = TRUE,
    quiet = TRUE,
    n_cores = parallel::detectCores() - 1,
    ...) {
        
        if (perm) {
            blocks_to_use=blocks
            for (k in seq(length(blocks)))
                blocks_to_use[[k]] <- as.matrix(blocks[[k]][sample(seq(NROW(blocks[[k]]))), ])
                rownames(blocks_to_use[[k]])=rownames(blocks[[k]])
        }else
        {
            blocks_to_use=blocks
        }
        args <- list(
            blocks = blocks_to_use,
            type = type,
            tol = tol,
            quiet = quiet,
            method = "complete",
            ...
        )
        
        args[[par]] <- par_value
        crit <- do.call(rgcca, args)$crit
        return(sum(sapply(crit, sum)))

}
