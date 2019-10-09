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

blocks.df <- lapply(keepCommonRow(blocks2[-5]), as.data.frame)
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

# Position of the NA values in clinic block
listNA <- which(is.na(blocks[[1]]), arr.ind = TRUE)
 
# Insert NA in the clinic post-processed
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
#       SGCCA
#####################

HV %in% sort(gsub("[a-zA-Z.]", "",  rownames(rna_dat)))

# Scaling
blocks.scaled <- scaling(blocks.df, TRUE)
load("blocks.scaled.RData")

sgcca.res <- rgcca.analyze(blocks = blocks.scaled, tau = c(1, 0.0735, 0.043, 0.143), scale = FALSE, type = "sgcca")
selected_var <- lapply(sgcca.res$a, function(x) which(x[, 1] != 0 | x[, 2] != 0))
sapply(selected_var, length)

lapply(selected_var, names)[-1]

updrs <- setResponse(blocks.scaled, "/home/etienne.camenen/DATA/Nucleiparks/UPDRS.tsv")
clusters <- setResponse(blocks.scaled, "/home/etienne.camenen/bin/fingerprint_clustering/clusters.tsv")
groups <- setResponse(blocks.scaled, "/home/etienne.camenen/DATA/Nucleiparks/group.tsv")

superblock <- Reduce(cbind, lapply(1:4, function(x) blocks.scaled[[x]][, selected_var[[x]]]))
write.table( imputeMean(superblock), file = "super_sgcca.tsv", sep="\t")

plotAVE(sgcca.res)

plotSamplesSpace(
                sgcca.res,
                updrs,
                1,
                2,
                1,
                TRUE,
                1,
                "UPDRS"
            )

plotSamplesSpace(
                sgcca.res,
                groups,
                1,
                2,
                1,
                TRUE,
                1,
                "Groups"
            )

plotSamplesSpace(
                sgcca.res,
                clusters,
                1,
                2,
                1,
                TRUE,
                1,
                "Clusters"
            )

plotVariablesSpace(
                rgcca = sgcca.res,
                blocks = list(superblock),
                comp_x = 1,
                comp_y = 2,
                superblock = TRUE,
                i_block = 1,
                text = TRUE,
                n_mark = 500,
                no_Overlap = TRUE,
                collapse = TRUE
            )

plotFingerprint(
        rgcca = sgcca.res, 
        blocks = blocks.scaled, 
        comp = 1, 
        #superblock = TRUE,
        type = "cor",
        n_mark = 500,
        i_block = 2,
        #collapse = TRUE
    )

sapply(blocks.scaled[-1], function(x) 1/sqrt(ncol(x)))

boot <- bootstrap(blocks = blocks.scaled, rgcca = sgcca.res, scale = FALSE, n_boot = 20000)
#boot2 <- boot[sample(1:length(boot), 500)]
selected_var <- getBootstrap(rgcca = sgcca.res, W = boot2, i_block = 2, collapse = FALSE)
#plotBootstrap(selected_var, sgcca.res, superblock = FALSE)

ggplot(selected_var, aes(x = abs(mean), y = weights, label = row.names(selected_var))) + geom_text(size=2) + labs(x = "Weights", y = "Non-zero occurences")

leaveOneOut.gcca(sgcca.res, A = blocks.scaled)
