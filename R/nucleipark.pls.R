library(RGCCA)
library(MASS)
library(ggplot2)
library(ggrepel)
library(visNetwork)

wd1 <- "../rgccaLauncher/R/"
for (f in c(list.files(wd1)[-1]))
    source(paste0(wd1, f))

funcs <- c("sgcca.crit", "sgcca.permute.crit", "rgcca", "rgccak", "sgcca", "sgccak", "defl.select", "initsvd", "pm", "scale3", "cov3", "norm2")
for (f in funcs)
    source(paste0("R/", f, ".R"))

#################
# Loading
#################

# files <- paste0(
#     "~/DATA/Nucleiparks/Nucleiparks_full/",
#     c("clinic", "metabolomic", "transcriptomic", "lipidomic") ,
#     ".tsv", 
#     collapse = ","
# )
# 
# blocks <- setBlocks(paste0(files, ","," /home/etienne.camenen/DATA/Nucleiparks/Nucleiparks_selectedVar/Imagery3.tsv"))
# names(blocks)[[5]] <- "imagery"
load(file = "blocks_nointer.RData")
blocks2 <- blocks

#################
# Cofounding
#################

omic = "imagery"
group <- setResponse(blocks.scaled, "/home/etienne.camenen/DATA/Nucleiparks/UPDRS.tsv")

blocks.df <- lapply(keepCommonRow(list(blocks2[["clinic"]], blocks2[[omic]])), as.data.frame)
# blocks.df <- lapply(blocks2, as.data.frame)
# blocks.df = list(blocks.df[["clinic"]])
blocks <- blocks.df
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

# # Position of the NA values in clinic block
listNA <- which(is.na(blocks[[1]]), arr.ind = TRUE)
# 
# # Insert NA in the clinic post-processed
for (i in listNA[,2]){
    blocks.df[[1]][[i]] <- append(
        blocks.df[[1]][[i]],
        NA,
        listNA[which(listNA[, 2] == i), 1] - 1
    )
}

# Convert in matrix
blocks.df[[1]] <- matrix(
    unlist(blocks.df[[1]]),
    nrow = nrow(blocks[[1]]),
    ncol = ncol(blocks[[1]]),
    dimnames = list(row.names(blocks[[1]]), colnames(blocks[[1]]))
)

# Remove cofounding variables
blocks.df[[1]] <- blocks.df[[1]][, -c(1:2, 4:6)]


#####################
#       PLS
#####################

# Scaling
blocks.scaled <- scaling(blocks.df, TRUE)
names(blocks.scaled) <- c("clinic", omic)

pls.res <- rgcca.analyze(blocks = blocks.scaled, scheme = "horst", tau = c(1, 1), scale = FALSE, type = "rgcca")

plotAVE(pls.res)

plotSamplesSpace(
                rgcca = pls.res,
                resp = group,
                comp_x = 1,
                comp_y = 2,
                i_block = 1,
                text = TRUE,
                reponse_name = "UPDRS",
                no_Overlap = FALSE
            )

superblock <- cbind(blocks.scaled[[1]], blocks.scaled[[2]][, which(pls.res$a[[2]][, 1] != 0 | pls.res$a[[2]][, 2] != 0)])
write.table( imputeMean(superblock), file = paste0("super_", omic, "_pls.tsv"), sep="\t")

group2 <- setResponse(blocks.scaled, "/home/etienne.camenen/bin/fingerprint_clustering/clusters.tsv")

plotSamplesSpace(
                rgcca = pls.res,
                resp = group2,
                comp_x = 1,
                comp_y = 2,
                i_block = 1,
                text = TRUE,
                reponse_name = "Clusters",
                no_Overlap = FALSE
            )

plotVariablesSpace(
                rgcca = pls.res,
                blocks = list(superblock),
                comp_x = 1,
                comp_y = 2,
                superblock = TRUE,
                i_block = 1,
                text = TRUE,
                n_mark = 50,
                collapse = TRUE,
                no_Overlap = TRUE
            )

plotFingerprint(
        rgcca = pls.res,
        blocks = blocks.scaled,
        comp = 1,
        i_block = 1,
        superblock = TRUE,
        type = "cor",
        n_mark = 30,
        collapse = TRUE
    )


#####################
#       PCA
#####################


blocks.super <- list(blocks.scaled[[1]], block2.selected, superblock)
names(blocks.super) <- c(names(blocks.scaled), "superblock")

J <- length(blocks.super)
connection <- matrix(0, J, J)
connection[seq_len(J - 1), J] <- connection[J, seq_len(J - 1)] <- 1

pls.super.res <- rgcca.analyze(blocks = blocks.super, connection = connection, scheme = "horst", scale = FALSE)

###
superblock <- imputeMean(pls.super.res$Y[["superblock"]])
# pca.super.res <- rgcca.analyze(blocks = list(superblock, superblock), scale = FALSE, type = "pca")

plotSamplesSpace(
                rgcca = pls.super.res,
                resp = group,
                comp_x = 1,
                comp_y = 2,
                i_block = 2,
                text = TRUE,
                reponse_name = "UPDRS"
            )

# J <- length(blocks.scaled)
# superblock <- c(as.list(1:J), list(Reduce(cbind,blocks.scaled)))
# pca.super.res$Y <- c(as.list(1:J), list(pca.super.res$Y[[1]]))
# pca.super.res$a <- c(lapply(blocks.scaled, colnames), list())
# vars <- colnames(superblock[[J + 1]])
# pca.super.res$a[[ J + 1]] <- matrix(Reduce(rbind, pls.res$a), length(vars), 2, dimnames = list(vars, 1:2))
# class(pca.super.res) <- "sgcca"

nodes <- getNodes(blocks.scaled, rgcca = pls.res)
edges <- getEdges(connection, blocks.scaled)
plotNetwork2(nodes, edges, blocks.scaled)



# group2 <- setResponse(blocks.scaled, "/home/etienne.camenen/bin/fingerprint_clustering/clusters.tsv")


plotVariablesSpace(
                rgcca = pls.super.res,
                blocks = blocks.super,
                comp_x = 1,
                comp_y = 2,
                superblock = TRUE,
                text = TRUE,
                n_mark = 100,
                no_Overlap = FALSE
            )

plotFingerprint(
        rgcca = pls.super.res, 
        blocks = blocks.super,
        comp = 1,
        superblock = FALSE,
        i_block = 2,
        type = "weight",
        n_mark = 100
    )

df <- getVar(pls.super.res, blocks.super, 1, 1, 2, "cor")
markers2 <- rownames(data.frame(getRankedValues(df, 1, TRUE), order = nrow(df):1)[1:30, ])
