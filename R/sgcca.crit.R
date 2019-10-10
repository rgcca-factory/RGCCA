# An intern function used by sgcca.permute to perform multiple sgcca with permuted rows

sgcca.crit <- function(A, C, c1s, ncomp, scheme, tol = .Machine$double.eps, scale = TRUE, perm = TRUE) {

    if(perm){
        for (k in 1:length(A))
            A[[k]] <- A[[k]][sample(1:nrow(A[[k]])),]
    }

    tasks <- simplify2array(parallel::mclapply(1:nrow(c1s),
        function(i) {
            out <- sgcca(
                    A = A,
                    C = C,
                    c1 = c1s[i,],
                    ncomp = ncomp,
                    scheme = scheme,
                    tol = tol,
                    scale = scale
                )
            return(c(mean(unlist(lapply(out$crit, function(x) x[length(x)]))), i))
        },
        mc.cores = 4))
    
    crit = rep(NA, NROW(c1s))
    
    for (i in seq(nrow(c1s)))
        crit[tasks[2, i]] <- tasks[1, i]

    return(crit)
}
