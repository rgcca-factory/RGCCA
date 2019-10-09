library(RGCCA)
library(MASS)

for (f in c("sgcca.crit", "sgcca", "sgccak", "defl.select", "initsvd", "pm", "scale3", "cov3", "norm2"))
    source(paste0("R/", f, ".R"))

parser <- argparse::ArgumentParser()
parser$add_argument("--job-id", required=T, metavar="integer")
parser$add_argument("--omic", required=T, metavar="character")
args <- parser$parse_args()

load(paste0(args$omic, ".RData"))

c1s <- seq(1 / sqrt(ncol(blocks.scaled[[1]])), 1, by = 0.01)
c1s <- as.matrix(cbind(c1s, c1s))

J <- length(blocks.scaled)
C <- matrix(0, J, J)
C[2:J, 1] <- C[1, 2:J] <- 1
    
sgcca.perm <- sgcca.crit(
    A = blocks.scaled,
    C = C,
    c1s = c1s,
    ncomp = rep(2, 2),
    scheme = "horst",
    scale = FALSE
)

write.table(
    as.matrix(t(sgcca.perm)),
    file = paste0(args$job_id, ".tsv"),
    col.names = FALSE,
    row.names = FALSE,
    sep = "\t"
)
