check_response <- function(response = NULL, df = NULL) {

    if (!is.null(response)) {
        qualitative <- is.character2(response)

        # if (length(qualitative) > 1)
        #     stop(
        #     "Please, select a response file with either qualitative data only or quantitative data only.",
        #     108
        #     )

        if (!qualitative)
            response <- to_numeric(response)

        if (NCOL(response) > 1) {
            disjunctive <- unique(apply(response, 1, sum))

            if (length(disjunctive) &&
                unique(disjunctive %in% c(0, 1)) && disjunctive) {

                response2 <- factor(apply(response, 1, which.max))

                if (!is.null(colnames(response)))
                    levels(response2) <- colnames(response)

                return(
                    as.matrix(
                        data.frame(
                            as.character(response2),
                            row.names = rownames(response)
                )))

            } else {
                warning("There is multiple columns in the response file. By default, only the first one is taken in account.")
                return(as.matrix(response[, 1]))
            }
        }

        return(response)

    } else
        return(rep(1, NROW(df[[1]])))

}
