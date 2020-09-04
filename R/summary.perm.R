summary.perm <- function(x) {
    tab_res <- cbind(x$crit, x$means, x$sds, x$zstat, x$pvals)
    colnames(tab_res) <- c("Crit", "Crit Perm", "Sd", "Z", "P-Value")
    round(cbind(x$penalties, tab_res), 3)
}
