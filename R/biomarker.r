biomarker <- function(
    resRGCCA,
    listData = NULL,
    block = 1,
    axes = 1,
    selec = "all",
    percentBm = 0.5) {

  
   A <- as.matrix(resRGCCA$call$blocks[[block]])
  
    
    x <- rep(NA, NCOL(A))

    for (i in 1:NCOL(A))
        x[i] <- cor(
                A[rownames(resRGCCA$Y[[block]]), i],
                as.matrix(resRGCCA$Y[[block]][rownames(resRGCCA$Y[[block]]), axes]),
                use = "pairwise.complete.obs"
            )

    names(x) <- colnames(A)
    if (selec == "all")
        selectionX <- 1:length(x)

    if (is.numeric(selec))
        selectionX <- 1:(min(selec, round(percentBm * length(x))))

    indices1 <- sort(
            abs(x),
            index.return = T,
            decreasing = TRUE,
            na.last = TRUE
        )$ix

    biomarkers <- rev((x[indices1])[selectionX])

    return(biomarkers)
}
