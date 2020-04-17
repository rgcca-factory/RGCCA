# Print and save indidvidual analysis attributes
# TODO: not only two components, "" for missing values
save_ind <- function(
    rgcca,
    compx = 1,
    compy = 2,
    file = "individuals.tsv") {
    
    stopifnot(is(rgcca, "rgcca"))
    for (i in c("compx", "compy")) {
        for (j in seq(length(rgcca$call$ncomp)))
            check_compx(i, get(i), rgcca$call$ncomp, j)
    }
    file <- paste0(file, collapse = " ")
    
    inds <- Reduce(cbind, lapply(
        rgcca$Y,
        function(x) x[, c(compx, compy)]))
    colnames(inds) <- as.vector(sapply(
        names(rgcca$call$blocks),
        function(x) paste0(x, ".axis", c(compx, compy))))

    write.table(inds, file, sep = "\t")

    invisible(inds)
}
