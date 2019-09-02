# An intern function used by sgcca.permute to perform multiple sgcca with permuted rows

sgcca.crit <- function(A, C, c1s, ncomp, scheme, tol, crit = crit, scale = TRUE) {

    for (k in 1:length(A))
        A[[k]] <- A[[k]][sample(1:nrow(A[[k]])), ]

    for (i in 1:NROW(c1s)) {
        out <- sgcca(A = A, C = C, c1 = c1s[i, ], ncomp = ncomp, scheme = scheme, tol = tol, scale = scale)
        crit[i] <-  mean(unlist(lapply(out$crit, function(x) x[length(x)])))
    }

    return(crit)
}
