source("load_nucleipark.R")

#################
# Best sparsity
#################

nperm <- 500
step <- .01
res.spls <- list()

for (i in 2:length(blocks)) {
    
    c1 <- seq(1/sqrt(ncol(blocks[[i]])), 1, by = step)
    c1s <- cbind(rep(1, length(c1)), c1)
    
    twoblocks <- scaling(list(blocks[[1]], blocks[[i]]), TRUE)
    
    res <- sgcca.permute.crit(
        A = twoblocks,
        c1s = c1s,
        nperm = nperm,
        scheme = "horst",
        ncomp = c(2,2)
    )

    res[["selected.variables"]] <- unique(names(which(res$sgcca.best$a[[2]][, 1] != 0 | res$sgcca.best$a[[2]][, 2] != 0)))
    plotFingerprint(res$sgcca.best, twoblocks, 1, F)

    res.spls[[i]] <- res
}

save(res.spls, file = "res3.RData")
