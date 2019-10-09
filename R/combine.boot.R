load("blocks.scaled.imagery.RData")

library(RGCCA)
library(MASS)

for (f in c("sgcca.crit", "sgcca", "sgccak", "defl.select", "initsvd", "pm", "scale3", "cov3", "norm2"))
    source(paste0("R/", f, ".R"))

c1s <- expand.grid(
    lapply(
        2:(length(blocks.scaled)-1), 
        function(x) seq(1 / sqrt(ncol(blocks.scaled[[x]])), 1, by = 0.1)
    )
)

c1s.1 <- rep(1, nrow(c1s))
c1s <- as.matrix(cbind(c1s.1, c1s, c1s.1))

J <- length(blocks.scaled)
C <- matrix(0, J, J)
C[2:J, 1] <- C[1, 2:J] <- 1

##############################

crits <- sgcca.crit(
    A = blocks.scaled,
    C = C,
    c1s = c1s,
    ncomp = rep(2, length(blocks.scaled)),
    scheme = "factorial",
    scale = FALSE,
    perm = FALSE
)

wd <- "/network/lustre/dtlake01/bioinfo/biostat/nucleipark/RGCCA/out_with_imaging/"
files <- c(list.files(wd))
files <- files[-length(files)]
permcrit <- c()

for (f in files)
    permcrit <- c(permcrit, read.table(paste0(wd, f), sep="\t"))

permcrit <- matrix(permcrit, nrow(c1s), length(files))

##############################

pvals <- zs <- matrix(NA, nrow = NROW(c1s), ncol = NCOL(c1s) + 1)


for (i in 1:NROW(c1s)) {
    pvals[i, ] <- c(c1s[i, ], mean(unlist(permcrit[i, ]) >= crits[i]))
    zs[i, ] <- c(c1s[i, ], (crits[i] - mean(unlist(permcrit[i, ]))) / (sd(unlist(permcrit[i, ]))))
}

bestpenalties <- c1s[which.max(zs[, length(blocks.scaled) + 1]), 1:length(blocks.scaled)]

sgcca.best <- sgcca(
        A = blocks.scaled,
        C = C,
        c1 = bestpenalties,
        ncomp = rep(2, length(blocks.scaled)),
        scheme = "factorial"
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

save(sgcca.res, file = "sgcca.perm.image.RData")
