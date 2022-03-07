# Print and save indidvidual analysis attributes
# TODO: not only two components, "" for missing values
save_ind <- function(rgcca,
                     file = "individuals.tsv") {
  stopifnot(is(rgcca, "rgcca"))
  file <- paste0(file, collapse = " ")

  inds <- inds <- Reduce(cbind, rgcca$Y)
  colnames(inds) <- unlist(sapply(
    names(rgcca$call$blocks),
    function(x) paste0(x, ".comp", seq(NCOL(rgcca$Y[[x]])))
  ))

  write.table(inds, file, sep = "\t")

  invisible(inds)
}
