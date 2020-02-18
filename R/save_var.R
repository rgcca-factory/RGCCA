# Print and save variables analysis attributes
save_var <- function(
    rgcca,
    compx = 1,
    compy = 2,
    file = "variables.tsv") {

    stopifnot(is(rgcca, "rgcca"))
    for (i in c("compx", "compy")) {
        for (j in seq(length(rgcca$call$ncomp)))
            check_compx(i, get(i), rgcca$call$ncomp, j)
    }
    file <- paste0(file, collapse = " ")

    indexes <- c("cor", "weight")

    vars <- Reduce(rbind, lapply(seq(length(rgcca$call$blocks)), function(i)
            data.frame(
                Reduce(cbind,
                        lapply(indexes, function(x)
                            get_ctr(rgcca, compx, compy, i_block = i, type = x))),
                names(rgcca$call$blocks)[i]
            )))

    colnames(vars) <- c(as.vector(sapply(indexes, function(x)
            paste0(x, ".", paste0("axis.", c(compx, compy))))), "block")

    write.table(vars, file, sep = "\t")

    invisible(vars)
}
