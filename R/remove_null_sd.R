#' Remove column having a standard deviation equals to 0
#'
#' @param list_m A list of dataframe
#' @return A list of dataframe
#' @examples
#' df = sapply(seq(3), function(x) runif(10))
#' df = cbind(df, rep(1, 10))
#' remove_null_sd(list(df))
#' @export
remove_null_sd <- function(list_m) {
    
    names <- names(list_m)
    
    column_sd_null <- lapply(list_m, 
                             function(x)
                                 which(apply(x, 2, sd) == 0))
    blocks_index <- seq(1, length(list_m))[
        unlist(
            lapply(
                column_sd_null,
                function(x) length(x) > 0))]
    
    list_m <- lapply(
        seq(length(list_m)),
        function(x) {
            if (x %in% blocks_index)
                list_m[[x]][, -column_sd_null[[x]]]
            else
                list_m[[x]]
        })
    
    names(list_m) <- names
    return(list_m)
}