#original data in rgcca
 # raw is the original data
 # blocks is the centered, scaled data with additional NA if a subject is missing
 # A is the imputed "blocks"
data(Russett)
blocks = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
              politic = Russett[, 6:11] )

resRgcca=rgcca(blocks)
resRgcca$call$raw[[1]][1,]
resRgcca$call$blocks[[1]][1,]
resRgcca$A[[1]][1,]

blocks_with_na=blocks
blocks_with_na[[1]][1,]=NA
blocks_with_na[[1]][2,2]=NA
resRgcca=rgcca(blocks_with_na,method = "mean")
resRgcca$call$raw[[1]][1,]
resRgcca$call$blocks[[1]][1,]
resRgcca$A[[1]][1,]

resRgcca=rgcca(blocks,superblock = T)
resRgcca$call$raw[[1]][1,]
resRgcca$call$blocks[[1]][1,]
resRgcca$A[[1]][1,]

lapply(resRgcca$call$raw,dim)
lapply(resRgcca$call$blocks,dim)
lapply(resRgcca$A,dim)

blocks_with_na=blocks
blocks_with_na[[1]][1,]=NA
blocks_with_na[[1]][2,2]=NA
resRgcca=rgcca(blocks_with_na,method = "mean",superblock=T)
resRgcca$call$raw[[1]][1,]
resRgcca$call$blocks[[1]][1,]
resRgcca$A[[1]][1,]
