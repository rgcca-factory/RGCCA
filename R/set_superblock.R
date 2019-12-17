set_superblock <- function(blocks, superblock = FALSE, type = "rgcca") {

    if (superblock | tolower(type) == "pca") {
        # if(type != 'pca') warnconnection('superblock')
        blocks[["superblock"]] <- Reduce(cbind, blocks)
    }
    return(blocks)
}
