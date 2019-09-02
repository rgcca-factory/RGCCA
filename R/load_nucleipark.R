library(RGCCA)
library(MASS)

wd1 <- "~/bin/rgccaLauncher/R/"
for (f in c(list.files(wd1)[-1]))
    source(paste0(wd1, f))

funcs <- c("sgcca.crit", "sgcca.permute.crit", "rgcca", "rgccak", "sgcca", "sgccak", "defl.select", "initsvd", "pm", "scale3", "cov3", "norm2")
for (f in funcs)
    source(paste0("~/bin/RGCCA/R/", f, ".R"))

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
