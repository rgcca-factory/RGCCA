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
    tau = rep(1, length(blocks)),
    sparsity = rep(1, length(blocks)),
    perm = TRUE,
    quiet = TRUE,
    n_cores = parallel::detectCores() - 1,
    ...) {

    if (perm) {
        for (k in seq(length(blocks)))
        blocks[[k]] <- as.matrix(blocks[[k]][sample(seq(NROW(blocks[[k]]))), ])
    }

    gcca <- quote(
    function(i) {
        func <- quote(
            rgcca(
                blocks = blocks,
                type = type,
                tol = tol,
                quiet = quiet,
                method = "complete",
                ...
            ))

        if(par[[1]] == "ncomp") {
            ncomp <- par[[2]][i, ]
            if (tolower(type) %in% c("sgcca", "spca", "spls"))
            func[["penalty"]] <- sparsity
            else
            func[["tau"]] <- tau
        } else
        func[[par[[1]]]] <- par[[2]][i, ]

        crit <- eval(as.call(func))$crit

        return(sum(sapply(crit, function(x) sum(x))))
    })

    varlist <- c()
    # get the parameter dot-dot-dot
    args_values <- c(...)
    # get the names of the arguments of function expect the ...
    args_func_names <- names(as.list(args("rgcca_crossvalidation")))
    # get only the names of the ... args
    args_dot_names <- setdiff(names(as.list(match.call()[-1])), args_func_names)
    n <- args_values
    if (!is.null(n))
        n <- seq(length(args_values))
    for (i in n) {
        # dynamically asssign these values
        assign(args_dot_names[i], args_values[i])
        # send them to the clusters to parallelize
        varlist <- c(varlist, args_dot_names[i])
        # without this procedure rgcca_crossvalidation(rgcca_res, blocks = blocks2)
        # or rgcca_crossvalidation(rgcca_res, blocks = lapply(blocks, scale)
        # does not work.
    }

    parallelize(
        varlist,
        seq(NROW(par[[2]])),
        eval(as.call(gcca)),
        n_cores = n_cores,
        envir = environment()
    )
}
