load("/network/lustre/dtlake01/bioinfo/biostat/nucleipark/data/blocks_without_cofounding/blocks.metlip.RData")

library(RGCCA)
library(MASS)

wd1 <- "../rgccaLauncher/R/"
for (f in c(list.files(wd1)[-1]))
    source(paste0(wd1, f))

for (f in c("crossvalidation.gcca", "predict.gcca", "sgcca", "sgccak", "defl.select", "initsvd", "pm", "scale3", "cov3", "norm2"))
    source(paste0("R/", f, ".R"))

parser <- argparse::ArgumentParser()
parser$add_argument("--job-id", required=T, metavar="integer")
args <- parser$parse_args()

step <- .01
c1s <- expand.grid(
    lapply(
        2:length(blocks.scaled),
        function(x) seq(1 / sqrt(ncol(blocks.scaled[[x]])) + step, 1, by = step)
    )
)

args = list(job_id=1)

c1s.1 <- rep(1, nrow(c1s))
c1s <- as.matrix(cbind(c1s.1, c1s))

J <- length(blocks.scaled)
C <- matrix(0, J, J)
C[2:J, 1] <- C[1, 2:J] <- 1

blocks.scaled[[1]] <- imputeMean(blocks.scaled[[1]])
blocks.scaled <- removeColumnSdNull(blocks.scaled)

sgcca.res <- sgcca(
    A = blocks.scaled,
    C = C,
    c1 = c1s[as.integer(args$job_id), ],
    ncomp = rep(2, J),
    scheme = "factorial",
    scale = FALSE
)
names(sgcca.res$a) <- names(blocks.scaled)

sgcca.cv <- crossvalidation.gcca(
    sgcca.res,
    A = blocks.scaled,
    bloc_to_pred = "clinic",
    validation = "loo"
)[[1]]

write.table(
    as.matrix(t(sgcca.cv)),
    file = paste0(args$job_id, ".tsv"),
    col.names = FALSE,
    row.names = FALSE,
    sep = "\t"
)
