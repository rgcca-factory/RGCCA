crits <- as.matrix(read.table("sgcca.tsv", sep = "\t"))
permcrit <- t(as.matrix(read.table("sgcca.permut.tsv", sep = "\t")))

c1s <- expand.grid(lapply(2:length(blocks.scaled), function(x) seq(1/sqrt(ncol(blocks.scaled[[x]])), 1, by = 0.05)))
c1s <- as.matrix(cbind(rep(1, nrow(c1s)), c1s))

pvals <- zs <- matrix(NA, nrow = NROW(c1s), ncol = NCOL(c1s) + 1)

# chaque ligne : [critère 1, critère 2, valeur pvals ou zs]
for (i in 1:NROW(c1s)) {
    pvals[i,] <- c(c1s[i,], mean(permcrit[i,] >= crits[i]))
    zs[i,] <- c(c1s[i,], (crits[i] - mean(permcrit[i,])) / (sd(permcrit[i,])))  # + 0.05))
}

bestpenalties <- c1s[which.max(zs[, length(A) + 1]), 1:length(A)]

sgcca.best <- sgcca(
        A,
        C = C,
        c1 = bestpenalties,
        ncomp = ncomp,
        scheme = scheme
    )

sgcca.res <- list(
        pvals = pvals,
        zstat = zs,
        bestpenalties = bestpenalties,
        sgcca.best = sgcca.best,
        permcrit = permcrit,
        crit = crits,
        penalties = c1s
    )
