print.perm <- function(x) {
    tab_res <- cbind(x$crit, x$means, x$sds, x$zstat, x$pvals)
    colnames(tab_res) <- c("RGCCA crit", "Perm. crit", "S.D.", "Z", "P-value")
    round(cbind(x$penalties, tab_res), 3)
}
