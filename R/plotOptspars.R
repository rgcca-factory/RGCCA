library(ggplot2)

omisc <- c( "metabolomic", "lipidomic", "transcriptomic")
load(paste0(omic, ".RData"))
load(paste0("pls.perm.", omic, ".RData"))

for (i in 1:length(omics)) {
    par.sv <- par()$mar
    par(mar = c(5.1, 5.1, 5.1, 2.1))
    zstat <- sgcca.res$zstat
    plot(
        zstat[, 2],
        zstat[, 3],
        type = "o",
        xlab = "",
        ylab = "",
        # xlab = "C1",
        # ylab = "Z-score",
        pch = "+",
        main = names(blocks.scaled)[1],
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
    
    sgcca.res$bestpenalties

# lapply(res.spls, function(x) length(x$selected.variables))
# lapply(res.spls, function(x) x$bestpenalties)
# 
# res.spls[[1]][["selected.variables"]] <- colnames(blocks[[1]])
# 
# blocks.sparsed <- mapply(
#     function(x, y) x[, y$selected.variables ],
#     blocks, 
#     res.spls
# )
