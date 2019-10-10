source("load_nucleipark.R")
load("sgcca.RData")

boot <- bootstrap(
    blocks.scaled,
    connection = C,
    c1 = sgcca.res$bestpenalties,
    ncomp = rep(2, 4),
    scheme = "factorial",
    scale = FALSE
)

save(boot, file = "bootstrap.RData")
