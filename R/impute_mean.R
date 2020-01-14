#' Impute NA
#'
#' Impute non available by means
#'
#' @param df A matrix containing non availables data
#' @examples
#' df = cbind(runif(9), runif(9))
#' df = rbind(df, c(NA, NA))
#' impute_mean(df)
#' @return A matrix with imputed values
#' @export
impute_mean <- function(df){
    
    # catch quantitative

    df <- to_numeric(df)

    if (any(is.na(df))) {
        df <- matrix(
            unlist(
                lapply(seq(NCOL(df)),
                function(x)
                    unlist(
                        lapply(as.list(df[, x]),
                        function(y)
                            ifelse(is.na(y),
                                mean(unlist(df[, x]), na.rm = TRUE), y))))),
            NROW(df),
            NCOL(df),
            dimnames = list(row.names(df), colnames(df)))
    }
    return(df)
}
