# Test on the sign of the correlation
check_sign_comp <- function(rgcca, w){

    w1 <- rgcca$a

    for (k in seq(length(w))) {
        for (j in seq(NCOL(w[[k]]))) {
            res <- cor(w1[[k]][, j], w[[k]][, j])
            if (!is.na(res) && res  < 0)
                w[[k]][, j] <- -1 * w[[k]][, j]
        }
    }

    return(w)
}
