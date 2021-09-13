get_comp_all <- function(
    rgcca,
    newA = rgcca$call$blocks,
    type = "train",
    newbloc_y = .Machine$integer.max,
    pred = NULL) {

    names <- unlist(lapply(rgcca$Y[-newbloc_y], colnames))

    if (type ==  "train")
        y <- lapply(
            seq(length(rgcca$Y)),
            function(x) rgcca$Y[[x]][row.names(rgcca$call$blocks[[x]]), ])
    else
        y <- pred

    # matrix of Y for all selected blocks and all components
    res  <- as.data.frame(Reduce("cbind", y[-newbloc_y]))

    # vector of character giving the name of the block and the number of the component
    colnames(res) <- paste(
        unlist(
            mapply(
                function(name, times)  rep(name, times),
                names(newA)[-newbloc_y],
                rgcca$call$ncomp
        )), names,
        sep = "_")

    return(res)
}
