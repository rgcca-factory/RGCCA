# An intern function used by sgcca.permute to perform multiple sgcca with permuted rows
# data("Russett")
# blocks = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
#     politic = Russett[, 6:11] )
# rgcca_permutation_k(blocks)
rgcca_permutation_k <- function(
    blocks,
    par = list("ncomp", expand.grid(rep(list(seq(2)), length(blocks)))),
    tol = 1e-03,
    type = "rgcca",
    sparsity = rep(1, length(blocks)),
    perm = TRUE,
    quiet = TRUE,
    n_cores = parallel::detectCores() - 1,
    ...) {

    if (perm) {
        for (k in seq(length(blocks)))
        blocks[[k]] <- as.matrix(blocks[[k]][sample(seq(NROW(blocks[[k]]))), ])
    }

    gcca <- function(i, ...) {
        
        
        args <- list(
            blocks = blocks,
            type = type,
            tol = tol,
            quiet = quiet,
            method = "complete"
        )
        
        inargs <- list(...)
        #args[[names(inargs)]] <- inargs
        args[[par[[1]]]] <- par[[2]][i, ]
        
        crit <- do.call(rgcca, args)$crit
        
        return(sum(sapply(crit, sum)))
    }

    varlist <- c()
    # get the parameter dot-dot-dot
    args_values <- list(...)
    args_names <- names(args_values)
    n <- args_values
    if (!is.null(n))
        n <- seq(length(args_values))
    for (i in n) {
        if (!is.null(args_names[i])) {
            # print("#")
            # print(args_names[i])
            # print(args_values[[i]])
            # dynamically asssign these values
            assign(args_names[i], args_values[[i]])
            # send them to the clusters to parallelize
            varlist <- c(varlist, args_names[i])
            # without this procedure rgcca_crossvalidation(rgcca_res, blocks = blocks2)
            # or rgcca_crossvalidation(rgcca_res, blocks = lapply(blocks, scale)
            # does not work.
        }
    }

    for (v in varlist)
        print(get(v))
    parallelize(
        varlist,
        seq(NROW(par[[2]])),
        gcca,
        n_cores = n_cores,
        envir = environment()
    )
}
