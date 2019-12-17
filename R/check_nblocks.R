check_nblocks <- function(blocks, type) {
    if (tolower(type) == "pca") {
        msg <- "Only one block is"
        exit_code <- 110
    } else{
        msg <- "Two blocks are"
        exit_code <- 111
    }

    stop(
        paste0(
            length(blocks),
            " blocks used in the analysis. ",
            msg ,
            " required for a ",
            type,
            "."
        ),
        exit_code = exit_code
    )
}
