check_method <- function(method) {
    analysis <- c("rgcca", "cpca-w", "gcca", "hpca", "maxbet-b", "maxbet", 
        "maxdiff-b","maxdiff", "maxvar-a", "maxvar-b", "maxvar", "niles", 
        "r-maxvar", "rcon-pca", "ridge-gca", "sabscor", "ssqcor", "ssqcor", 
        "ssqcov-1", "ssqcov-2", "ssqcov", "sum-pca", "sumcor", "sumcov-1", 
        "sumcov-2", "sumcov", "sabscov", "plspm", "cca", "ra", "ifa", "pls",
        "pca", "sgcca", "spls", "spca")
    if (!tolower(method) %in% analysis)
        stop(
            paste0("Wrong type of analysis. Please select one among the following
            list: ", paste(analysis, collapse = ", ")),
            exit_code = 112
        )
}