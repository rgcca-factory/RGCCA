load("blocks.scaled.RData")

library(RGCCA)
library(MASS)

for (f in c("sgcca.crit", "sgcca", "sgccak", "defl.select", "initsvd", "pm", "scale3", "cov3", "norm2"))
    source(paste0("R/", f, ".R"))

c1s <- expand.grid(
    lapply(
        2:length(blocks.scaled), 
        function(x) seq(1 / sqrt(ncol(blocks.scaled[[x]])), 1, by = 0.1)
    )
)
c1s <- as.matrix(cbind(rep(1, nrow(c1s)), c1s))

J <- length(blocks.scaled)
C <- matrix(0, J, J)
C[2:J, 1] <- C[1, 2:J] <- 1
    
sgcca.perm <- sgcca.crit(
    A = blocks.scaled,
    C = C,
    c1s = c1s,
    ncomp = rep(2, 4),
    scheme = "factorial",
    scale = FALSE
)

write.table(
    as.matrix(t(sgcca.perm)),
    file = "sgcca.permut.tsv",
    append = TRUE,
    col.names = FALSE,
    row.names = FALSE,
    sep = "\t"
)
