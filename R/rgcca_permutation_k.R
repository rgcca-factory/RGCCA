# An intern function used by sgcca.permute to perform multiple sgcca with permuted rows
# data("Russett")
# blocks = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
#     politic = Russett[, 6:11] )
# rgcca_permutation_k(blocks)
rgcca_permutation_k <- function(
    blocks,
    connection=NULL,
    par = "ncomp",
    par_value=rep(1,length(blocks)),
    tol = 1e-03,
    type = "rgcca",
    sparsity = rep(1, length(blocks)),
    tau=rep(1,length(blocks)),
    perm = TRUE,
    quiet = TRUE,
    n_cores = parallel::detectCores() - 1,
    superblock=FALSE,
    scale=TRUE,
    scale_block=TRUE,
    scheme="factorial",
    method="nipals",
    ncomp=rep(1,length(blocks)),
    bias=FALSE) {
     t0=Sys.time()
        if (perm) {
            blocks_to_use=blocks
            blocks_to_use=lapply(seq(length(blocks)),function(k)
                { 
                    blocks_to_use_k <- as.matrix(blocks[[k]][sample(seq(NROW(blocks[[k]]))), ])
                    rownames(blocks_to_use_k)=rownames(blocks[[k]])
                    return(blocks_to_use_k)
                })
            names(blocks_to_use)=names(blocks)

            
          #  for (k in seq(length(blocks)))
          #      blocks_to_use[[k]] <- as.matrix(blocks[[k]][sample(seq(NROW(blocks[[k]]))), ])
          #      rownames(blocks_to_use[[k]])=rownames(blocks[[k]])
        }else
        {
            blocks_to_use=blocks
        }

        if(par=="ncomp")
        {
            res <- rgcca(
                blocks = blocks_to_use,
                type = type,
                tol = tol,
                quiet = quiet,
                method = method,
                superblock=superblock,
                scale=scale,
                scale_block=scale_block,
                scheme=scheme,
                connection=connection,
                ncomp=par_value,
                sparsity=sparsity,
                tau=tau
            )
        }
        if(par=="tau")
        {
            res <- rgcca(
                blocks = blocks_to_use,
                type = type,
                tol = tol,
                quiet = quiet,
                method = method,
                superblock=superblock,
                scale=scale,
                scale_block=scale_block,
                scheme=scheme,
                connection=connection,
                ncomp=ncomp,
                sparsity=NULL,
                tau=par_value
            )
        }
        if(par=="sparsity")
        {
            res <- rgcca(
                blocks = blocks_to_use,
                type = type,
                tol = tol,
                quiet = quiet,
                method = method,
                superblock=superblock,
                scale=scale,
                scale_block=scale_block,
                scheme=scheme,
                connection=connection,
                ncomp=ncomp,
                sparsity=par_value,
                tau=NULL
            )
        }
        crit <- res$crit[length(res$crit)]
        return(sum(crit))
}
