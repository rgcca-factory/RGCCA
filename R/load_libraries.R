load_libraries <- function(librairies) {
    for (l in librairies) {
        if (!(l %in% installed.packages()[, "Package"]))
            utils::install.packages(l, repos = "http://cran.wustl.edu")
        library(
            l,
            character.only = TRUE,
            warn.conflicts = FALSE,
            quietly = TRUE
        )
    }
}
