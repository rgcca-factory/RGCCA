source("load_nucleipark.R")

#################
# Best sparsity
#################

# nperm <- 500
nperm <- 2
step <- .5

c1s <- expand.grid(lapply(2:length(blocks.scaled), function(x) seq(1/sqrt(ncol(blocks.scaled[[x]])), 1, by = step)))
c1s <- as.matrix(cbind(rep(1, nrow(c1s)), c1s))

J <- length(blocks.scaled)
C <- matrix(0, J, J)
C[2:J, 1] <- C[1, 2:J] <- 1

    
sgcca.res <- sgcca.permute.crit(
    A = blocks.scaled,
    C = C,
    c1s = c1s,
    nperm = nperm,
    ncomp = rep(2, 4),
    scheme = "factorial",
    scale = FALSE
)

# sgcca.res[["selected.variables"]] <- lapply(2:unique(names(which(sgcca.res$sgcca.best$a[[2]][, 1] != 0 | sgcca.res$sgcca.best$a[[2]][, 2] != 0)))
# plotFingerprint(sgcca.res$sgcca.best, twoblocks, 1, F)

save(sgcca.res, file = "sgcca.RData")
