load(file = "sgcca_500_.02.RData")
library(ggplot2)

names(sgcca.res$sgcca.best$a) <- names(blocks.scaled)
plotFingerprint(sgcca.res$sgcca.best, blocks.scaled, superblock = TRUE, type = "weight", collapse = TRUE)

J <- length(blocks.scaled)
C <- matrix(0, J, J)
C[2:J, 1] <- C[1, 2:J] <- 1

boot <- bootstrap(
    blocks.scaled,
    n_boot = 14,
    connection = C,
    tau = sgcca.res$bestpenalties,
    ncomp = rep(2, 4),
    scheme = "factorial",
    scale = FALSE,
    type = "sgcca"
)

biomarkers <- getBootstrap(sgcca.res$sgcca.best, boot, comp = 1, collapse = TRUE)
plotBootstrap(biomarkers, sgcca.res$sgcca.best, TRUE, n_mark = 30)
