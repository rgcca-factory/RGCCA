# #' @export
set_superblock <- function(blocks, superblock = FALSE, method = "rgcca", verbose =  TRUE) {

    if (superblock | tolower(method) == "pca") {
        if (tolower(method) != 'pca' && verbose)
            warn_connection('superblock')
        blocks[["superblock"]] <- Reduce(cbind, blocks)
    }
    return(blocks)
}
