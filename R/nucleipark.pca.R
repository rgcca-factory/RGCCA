library(RGCCA)
library(MASS)
library(ggplot2)
library(ggrepel)

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
omic <- "clinic"
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

# Position of the NA values in clinic block
listNA <- which(is.na(blocks[[1]]), arr.ind = TRUE)
#
# # Insert NA in the clinic post-processed
for (i in listNA[,2])
    blocks.df[[1]][[i]] <- append(
        blocks.df[[1]][[i]],
        NA,
        listNA[which(listNA[, 2] == i), 1] - 1
    )

# # Convert in matrix
blocks.df[[1]] <- matrix(
    unlist(blocks.df[[1]]),
    nrow = nrow(blocks[[1]]),
    ncol = ncol(blocks[[1]]),
    dimnames = list(row.names(blocks[[1]]), colnames(blocks[[1]]))
)


# Remove cofounding variables
blocks.df[[1]] <- blocks.df[[1]][, -c(1:2, 4:6)]
blocks <- list(blocks.df[[2]])

# Scaling
blocks.scaled <- scaling(blocks, TRUE)
blocks.scaled = list(blocks.scaled[[1]], blocks.scaled[[1]])
names(blocks.scaled) <- c(omic, "superblock")
# save(blocks.scaled, file = "transcriptomic.RData")

load(paste0(omic, ".RData"))
load(paste0("pls.perm.", omic, ".RData"))

pca.res <- rgcca.analyze(blocks.scaled, scale = FALSE, type = "pca")
#pca.res <- sgcca.res$sgcca.best
names(pca.res$a) <- names(blocks.scaled)

selectedVar <- which(pca.res$a[[2]][, 1] != 0 | pca.res$a[[2]][, 2] != 0)
length(selectedVar)
write.table( imputeMean(blocks.scaled[[2]][, selectedVar]), file = paste0(omic, "_cofunding.tsv"), sep="\t")


plotAVE(pca.res)
updrs <- setResponse(blocks.scaled, "/home/etienne.camenen/DATA/Nucleiparks/UPDRS.tsv")
clusters <- setResponse(blocks.scaled, "/home/etienne.camenen/bin/fingerprint_clustering/clusters.tsv")
groups <- setResponse(blocks.scaled, "/home/etienne.camenen/DATA/Nucleiparks/group.tsv")
plotSamplesSpace(
                pca.res,
                updrs,
                1,
                2,
                1,
                TRUE,
                1,
                "UPDRS"
            )

plotSamplesSpace(
                pca.res,
                groups,
                1,
                2,
                1,
                TRUE,
                1,
                "Groups"
            )

plotSamplesSpace(
                pca.res,
                clusters,
                1,
                2,
                1,
                TRUE,
                1,
                "Clusters"
            )


corResponse(pca.res,
    blocks.df,
    i_response = 1,
    comp = 1,
    i_block = 1)

plotVariablesSpace(
                rgcca = pca.res,
                blocks = blocks.scaled,
                comp_x = 1,
                comp_y = 2,
                superblock = FALSE,
                text = TRUE,
                n_mark = 100,
                no_Overlap = TRUE
            )

    plotFingerprint(
        rgcca = pca.res, 
        blocks = blocks.scaled, 
        comp = 1, 
        superblock = FALSE,
        type = "cor",
        n_mark = 30,
    )




blocks3 <- blocks2[1:4]

library(VennDiagram)

venn <- function(blocks){
    x11()
    
    commonRows <- function(i_block_1, i_block_2) {
        if (length(i_block_2) < 2)
            i_block_2 <- row.names(blocks[[i_block_2]])
        intersect(row.names(blocks[[i_block_1]]), i_block_2)
    }
    
    J <- length(blocks)
    n_rows <- sapply(blocks, nrow)
    f <- c("pairwise", "triple", "quad", "quintuple")[J-1]
    
    if ( J < 2 || J > 5)
        stop("The number of blocks must be between 2 and 5.")
        
    func <- quote(
        get(paste0("draw.", f, ".venn"))(
            category = names(blocks[1:J]),
            fill = colorGroup(names(blocks)[1:J]),
            cat.col = colorGroup(names(blocks)[1:J])
        )
    )
    
    if ( J == 2 )
        func$cross.area = length(commonRows(1,2))
    
    for (i in 1:J)
        func[ paste0("area", i) ] <- n_rows[i]
    
    if (J > 2){
        comb2 <- combn(1:J, 2)
        for (i in 1:ncol(comb2))
            func[paste0('n', paste(comb2[,i], collapse = ''))] <- length(commonRows(comb2[1, i], comb2[2, i]))
        comb3 <- combn(1:J, 3)
        for (i in 1:ncol(comb3))
            func[paste0('n', paste(comb3[,i], collapse = ''))] <- length(commonRows(comb3[1, i], commonRows(comb3[2, i], comb3[3, i])))
    }
    if (J > 3){
        comb4 <- combn(1:J, 4)
        for (i in 1:ncol(comb4))
            func[paste0('n', paste(comb4[,i], collapse = ''))]  <- length(intersect(commonRows(comb4[1, i], comb4[2, i]), commonRows(comb4[3, i], comb4[4, i])))
    }
    if (J == 5 )
        func$n12345 <- length(commonRows(1, intersect(commonRows(2, 3), commonRows(4, 5))))
    
    invisible(eval(as.call(func)))
}
