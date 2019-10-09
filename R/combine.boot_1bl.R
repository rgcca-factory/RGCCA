library(RGCCA)
library(MASS)

for (f in c("sgcca.crit", "sgcca", "sgccak", "defl.select", "initsvd", "pm", "scale3", "cov3", "norm2"))
    source(paste0("R/", f, ".R"))


for (omic in c("metabolomic", "lipidomic", "transcriptomic")){
    
  load(paste0(omic, ".RData"))

    c1s <- seq(1 / sqrt(ncol(blocks.scaled[[1]])), 1, by = 0.01)
    c1s <- as.matrix(cbind(c1s, c1s))
    
    J <- length(blocks.scaled)
    C <- matrix(0, J, J)
    C[2:J, 1] <- C[1, 2:J] <- 1
    
    ##############################
    
    crits <- sgcca.crit(
        A = blocks.scaled,
        C = C,
        c1s = c1s,
        ncomp = rep(2, length(blocks.scaled)),
        scheme = "horst",
        scale = FALSE,
        perm = FALSE
    )
    
    print(omic)
        
    wd <- paste0("/network/lustre/dtlake01/bioinfo/biostat/nucleipark/RGCCA/out/", omic, "/")
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
            scheme = "horst",
            scale = FALSE
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
    
    save(sgcca.res, file = paste0("pls.perm.2.", omic, ".RData"))

}
