# Print and save variables analysis attributes
save_var <- function(
    rgcca,
    compx = 1,
    compy = 2,
    file = "variables.tsv") {

    indexes <- c("cor", "weight")

    vars <- Reduce(rbind, lapply(seq(length(rgcca$blocks)), function(i)
            data.frame(
                Reduce(cbind,
                        lapply(indexes, function(x)
                            get_ctr(rgcca, compx, compy, i_block = i, type = x))),
                names(rgcca$blocks)[i]
            )))

    colnames(vars) <- c(as.vector(sapply(indexes, function(x)
            paste0(x, ".", paste0("axis.", c(compx, compy))))), "block")

    write.table(vars, file, sep = "\t")

    invisible(vars)
}
