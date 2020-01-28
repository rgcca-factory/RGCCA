# #' @export
set_superblock <- function(blocks, superblock = FALSE, type = "rgcca", verbose =  TRUE) {

    if (superblock | tolower(type) == "pca") {
        if (tolower(type) != 'pca' && verbose)
            warnconnection('superblock')
        blocks[["superblock"]] <- Reduce(cbind, blocks)
    }
    return(blocks)
}
