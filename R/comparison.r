#'Compares 2 RGCCA
#'
#' @param rgcca1 A result of a RGCCA function 
#' @param rgcca2 Another result of RGCCA function
#' @param nAxe=1 Number of axes taken into account during the comparison
#' @param selec Number of biomarkers to be selected
#' @param selectPatient A vector allowing to select only some patients for the RV calculation
#' @return \item{A} A list containing: a: the correlation between axes, rv: the rv coefficient, bm the biomarkers
#' @return \item{crit} Convergence criterion : abs(1-obj_k/obj_{k-1})
#' @return \item{obj} Vector containing the mean square error between the predict values and the original non missing values at each iteration
#' @title comparison of two RGCCA results
#' @examples 
#'  data();...
#'  

comparison <- function(
    rgcca1,
    rgcca2,
    nAxe = 1,
    selec = 10,
    selectPatient = NULL) {


    diffNorm2 <- function(vec1, vec2) {
        if (!is.null(vec1) & !is.null(vec2)) {
            c1 <- cor(vec1, vec2)
            c2 <- cor(vec1, -vec2)
            return(max(c1, c2))
        } else {
            return(NA)
        }
    }

    if (is.null(selectPatient)) {
        selectPatient <- rownames(rgcca1[["Y"]][[1]])
    }

    J <- length(rgcca1$A)
    com <- rv <- pctBm <- rep(NA, J)

    for (i in 1:J) {
        rv[i] <- FactoMineR::coeffRV(
                rgcca1[["Y"]][[i]][selectPatient, 1:2], 
                rgcca2[["Y"]][[i]][selectPatient, 1:2]
            )$rv
        com[i] <- diffNorm2(
                rgcca1[["astar"]][[i]][, nAxe], 
                rgcca2[["astar"]][[i]][, nAxe]
            )
        bm <- lapply(list(rgcca1, rgcca2), 
            function(x){
                biomarker(
                    resRGCCA = x,
                    block = i,
                    axes = nAxe,
                    selec = selec
                )
            }
        )
        pctBm[i] <- sum(names(bm[[2]]) %in% names(bm[[1]])) / length(bm[[1]])
    }
    
    return(list(a = com, rv = rv, bm = pctBm))
}
