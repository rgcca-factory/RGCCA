load_libraries <- function(librairies) {
    for (l in librairies) {
        if (!(l %in% installed.packages()[, "Package"]))
            utils::install.packages(l, repos = "http://cran.us.r-project.org")
        library(
            l,
            character.only = TRUE,
            warn.conflicts = FALSE,
            quietly = TRUE
        )
    }
}
