# Print and save variables analysis attributes
save_var <- function(
    rgcca,
    file = "variables.tsv") {

    stopifnot(is(rgcca, "rgcca"))
    file <- paste0(file, collapse = " ")

    indexes <- c("loadings", "weight")

    vars <- lapply(
        seq(length(rgcca$call$blocks)),
        function(i) {
            df <- data.frame(
                Reduce(cbind,
                        lapply(indexes, function(x)
                            Reduce(cbind, lapply(seq(rgcca$call$ncomp[i]), function(j) get_ctr(rgcca, j, j, i_block = i, type = x)[, 1, drop = FALSE]))
                            )),
                names(rgcca$call$blocks)[i]
            )
            colnames(df) <- c(unlist(lapply(
                indexes,
                function(x) sapply(seq(rgcca$call$ncomp[[i]]), function(y) paste0(x, ".", y)))), "blocks")
            df
        }
    )

    common <- unique(unlist(lapply(vars, names)))
    vars <- Reduce(
        rbind,
        lapply(vars,
               function(x) {
                   x[, setdiff(common, names(x))] <- NA
                   return(x)
                }))

    write.table(vars, file, sep = "\t")

    invisible(vars)
}
