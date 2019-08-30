source("~/bin/rgccaLauncher/R/parsing.R")
source("~/bin/rgccaLauncher/R/select.type.R")
source("~/bin/rgccaLauncher/R/plot.R")

source("~/bin/RGCCA/R/sgcca.permute.crit.R")
source("~/bin/RGCCA/R/sgcca.crit.R")

# rgcca
# rgccak
# defl.select
# initsvd

library(RGCCA)
library(ggplot2)
library(MASS)

source("~/bin/RGCCA/R/pm.R")
source("~/bin/RGCCA/R/scale3.R")
source("~/bin/RGCCA/R/cov3.R")
source("~/bin/RGCCA/R/norm2.R")

#################
# Loading
#################

blocks <- setBlocks(
    paste0(
        "~/DATA/Nucleiparks/Nucleiparks_full/",
        c("clinic", "metabolomic", "transcriptomic", "lipidomic") ,
        ".tsv", 
        collapse = ","
        )
    )

#################
# Missing data
#################

source("~/bin/RGCCA/R/intersection.R")
source("~/bin/RGCCA/R/MIRGCCA.R")
source("~/bin/RGCCA/R/imputeNN.R")
source("~/bin/RGCCA/R/imputeSB.R")
source("~/bin/RGCCA/R/biomarker.r")
source("~/bin/RGCCA/R/plotAnalysis.r")

library(PMA)
library(nipals)
library(FactoMineR)
library(missMDA)

nbSimul <- 5

# TODO: crit[iter] = NA when na.rm = FALSE
refRgcca <- rgcca(
    A = blocks,
    ncomp = rep(2,length(blocks)),
    scheme = "factorial",
    verbose = FALSE,
    returnA = TRUE
)

lapply(block_simul, function(x) unique(as.vector(unique(apply(x, 1, is.na)))))

res <- parallel::mclapply(1:nbSimul, mc.cores = parallel::detectCores() - 1, FUN = function(x) {
    
    pNA <- 0.05
    listRgcca <- list()
    block_simul <- blocks
    
    for (i in 1:length(block_simul)) {

        n <- nrow(block_simul[[i]])
        p <- ncol(block_simul[[i]])
        nbNA <- round(pNA * n * p)
        listeIndicePossible <- merge(1:n, 1:p)
        indToRemove <- sample(1:(n * p), nbNA, replace = FALSE)
        listeIndicesToRemove <- listeIndicePossible[indToRemove, ]

        for (j in 1:length(indToRemove))
            block_simul[[i]][listeIndicesToRemove[j, "x"], listeIndicesToRemove[j, "y"]] <- NA
    }

    
    listRgcca[["MI-kNNAll"]] <- MIRGCCA(
        block_simul,
        k = "all",
        ni = 5,
        output = "weightedMean",
        scheme = "factorial",
        returnA = TRUE
    )$rgcca0

    block_simul_SB <- imputeSB(
        block_simul,
        ni = 20,
        ncomp = rep(2,length(blocks)),
    )
    
    listRgcca[["EM"]] <- rgcca(
        A = block_simul_SB$A,
        ncomp = rep(2,length(blocks)),
        scheme = "factorial",
        verbose = FALSE,
        returnA = TRUE
    )
    
    listRgcca[["Nipals"]] <- rgcca(
        block_simul,
        ncomp = rep(2,length(blocks)),
        scheme = "factorial",
        verbose = FALSE,
        returnA = TRUE,
        na.rm = TRUE
    )
    
    lapply(listRgcca, function(x)
        comparison(
            rgcca1 = refRgcca,
            rgcca2 = x
        )
    )

})

plotAnalysis(res, ylim = c(0, 1.1), output = "rv")
plotAnalysis(res, ylim = c(0, 1.1), output = "bm")
plotAnalysis(res, ylim = c(0, 1.1), output = "a")

block_simul_SB <- imputeSB(
    blocks,
    ni = 20,
    ncomp = rep(2,length(blocks)),
)


#################
# Cofounding
#################

blocks.df <- lapply(blocks, as.data.frame)
cl <- blocks.df[[1]]

# Weight by the cofunding effect residuals
blocks.df <- lapply(
    blocks.df, 
    function(x) apply(
            x,
            2, 
            function(y) 
                lm(y ~ cl$gender + cl$age + cl$weight + cl$height + cl$BMI, na.action = "na.exclude")$residuals
        )
)

# Position of the NA values in clinic block
listNA <- which(is.na(blocks[[1]]), arr.ind = TRUE)

# Insert NA in the clinic post-processed
for (i in listNA[,2])
    blocks.df[[1]][[i]] <- append(
        blocks.df[[1]][[i]],
        NA,
        listNA[which(listNA[, 2] == i), 1] - 1
    )

# Convert in matrix
blocks.df[[1]] <- matrix(
    unlist(blocks.df[[1]]),
    nrow = nrow(blocks[[1]]),
    ncol = ncol(blocks[[1]]),
    dimnames = list(row.names(blocks[[1]]), colnames(blocks[[1]]))
)


# Remove cofounding variables
blocks.df[[1]] <- blocks.df[[1]][, -c(1:2, 4:6)]
blocks <- blocks.df

#################
# Best sparsity
#################

nperm <- 500
step <- .01
res.spls <- list()

# expand.grid(lapply(2:length(blocks), function(x) seq(1/sqrt(ncol(blocks[[x]])), 1, by = .1)))

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

for (i in 2:length(blocks)) {
    #png(file = paste0(names(blocks)[i], "-bestc1.png"))
    par.sv <- par()$mar
    par(mar = c(5.1, 5.1, 5.1, 2.1))
    zstat <- res.spls[[i]]$zstat
    plot(
        zstat[, 2],
        zstat[, 3],
        type = "o",
        xlab = "",
        ylab = "",
        # xlab = "C1",
        # ylab = "Z-score",
        pch = "+",
        main = names(blocks)[i],
        axes = F,
        cex.lab = 1.5, font.lab = 3, font.axis = 3, cex.axis = 0.8, cex.main = 2, cex = 1, lwd = 3
    )
    points(
        zstat[which.max(zstat[, 3]), 2],
        zstat[which.max(zstat[, 3]), 3],
        cex = 3,
        col = "red",
        pch = "+"
    )
    abline(
        h = c(1.96, 2.58, 3.29),
        lty = 2,
        col = "grey",
        lwd = 4
    )
    axis(1, lwd = 4)
    axis(2, lwd = 4)
    par(par.sv)
    #dev.off()
}

lapply(res.spls, function(x) length(x$selected.variables))
lapply(res.spls, function(x) x$bestpenalties)

res.spls[[1]][["selected.variables"]] <- colnames(blocks[[1]])

blocks.sparsed <- mapply(
    function(x, y) x[, y$selected.variables ],
    blocks, 
    res.spls
)


#################
# RGCCA
#################

blocks.imputed <- c(clinic = list(imputeMean(blocks.sparsed[[1]])), blocks.sparsed[c(2:4)])
# t = lapply(blocks.imputed, dim)
# matrix(unlist(t), 2, 4, dimnames= list(c("r", "c"), colnames = names(t)))

blocks.intersected <- lapply(blocks.sparsed, function(x) x[ -c(unique(listNA[,1 ])), ])

# blocks

J <- length(blocks)
C <- matrix(0, J, J)
C[2:J, 1] <- C[1, 2:J] <- 1

best.tau <- rgcca(
   blocks.imputed,
    C = C,
    tau = "optimal",
    ncomp = rep(2, 4),
    scheme = "factorial", 
    verbose = F
)$tau

round(best.tau, 2)

##

blocks.sparsed.scaled <- scaling(blocks.sparsed, TRUE)

rgcca.res <- rgcca(
    blocks.sparsed.scaled,
    C = C,
    tau = best.tau,
    ncomp = rep(2, 4),
    scheme = "factorial",
    verbose = F,
    scale = FALSE
)

plotFingerprint(
    rgcca = rgcca.res,
    blocks = blocks.sparsed.scaled,
    superblock = FALSE,
    i_block = 4,
    n_mar = 30,
    type = "cor"
)

boot <- bootstrap(
    blocks.sparsed.scaled,
    connection = C,
    tau = best.tau,
    ncomp = rep(2, 4),
    scheme = "factorial",
    scale = FALSE
)

plotBootstrap(boot, n_mark = 30, i_block = 4)

# C = "superblock",
# matrix(c(best.tau, 0, 0), 2, 5)

#################
# PCA
#################

library(ade4)

pca <- dudi.pca(
    Reduce(cbind, rgcca.res$Y),
    scannf = FALSE,
    nf = 2
)

pca$AVE$AVE_X <- c(as.list(1:4), list(pca$eig/sum(pca$eig)))
pca$Y <- c(as.list(1:4), list(pca$li))
pca$a <- c(lapply(blocks.sparsed.scaled, colnames), list())

superblock <- c(as.list(1:4), list(Reduce(cbind,blocks.sparsed.scaled)))
vars <- colnames(superblock[[5]])
pca$a[[5]] <- matrix(c(vars, vars), length(vars), 2, dimnames = list(vars, 1:2))

fig <- plotFingerprint(rgcca = pca, blocks = superblock, comp= 2, n_mar = 30)
fig <- plotVariablesSpace(rgcca = pca, blocks = superblock)
fig <- plotSamplesSpace(rgcca = pca, resp = setResponse(file = "~/DATA/Nucleiparks/UPDRS.tsv"), reponse_name = "UPDRS")
savePlot("fig.png", fig)

pca <- rgcca(
    A = Reduce(cbind, rgcca.res$Y),
    scheme = "horst",
    ncomp = c(2,2)
    )
