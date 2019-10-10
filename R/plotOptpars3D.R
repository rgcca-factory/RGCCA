library(plotly)

load("sgcca.perm.whithout.image.Rdata")
zstat <- as.data.frame(sgcca.res$zstat[, -1])
colnames(zstat) <- c(paste0("block_", 2:4), "z")

sign <- as.double(zstat$z > qnorm(1 - 0.05 / 2)) + as.double(zstat$z > qnorm(1 - 0.01 / 2)) + as.double(zstat$z > qnorm(1 - 0.001 / 2)) 

plotly::plot_ly(
    zstat,
    x = ~ block_2,
    y = ~ block_3,
    z = ~ block_4,
    marker = list(
        color = ~ sign,
        showscale = TRUE,
        colorbar = list(
            title = 'Z-score'
        ),
        colorscale = list(
            list(0, "rgb(165,0,38)"), list(mean(zstat$z), "rgb(254,224,144)"), list(max(zstat$z), "rgb(49,54,149)")
        ),
        cauto = F,
        cmin = 0,
        cmax = max(zstat$z)
    )
) %>% add_markers()

zstat[which(zstat$z > qnorm(1 - 0.05 / 2)), ]


getRankedValues(zstat, comp = 4)

thresh <- 0.3
getRankedValues(zstat[which(zstat$z > 1.96 & zstat$block_2 < thresh & zstat$block_3 < thresh & zstat$block_4 < thresh), ], 4)

sgcca.res <- sgcca(
        A = blocks.scaled,
        C = C,
        c1 = c(1, 0.07, 0.041, 0.14),
        ncomp = rep(2, length(blocks.scaled)),
        scheme = "factorial"
    )

markers <- lapply(2:length(blocks.scaled), function(x) which(sgcca.res$a[[x]][, 1] != 0 | sgcca.res$a[[x]][, 2] != 0))
lapply(markers, length)
