#' Extract a bootstrap
#'
#' Extract statistical information from a bootstrap
#'
#' @inheritParams bootstrap
#' @inheritParams plot_histogram
#' @inheritParams plot_var_2D
#' @inheritParams plot_var_1D
#' @param w A list of list weights (one per bootstrap per blocks)
#' @return A matrix containing the means, 95% intervals, bootstrap ratio and p-values
#' @examples
#' library(RGCCA)
#' data("Russett")
#' blocks = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
#'     politic = Russett[, 6:11] )
#' rgcca_out = rgcca.analyze(blocks)
#' boot = bootstrap(rgcca_out, 2, FALSE, n_cores = 1)
#' get_bootstrap(rgcca_out, boot, n_cores = 1)
#' @export
get_bootstrap <- function(
    rgcca,
    w,
    comp = 1,
    i_block = length(w[[1]]),
    collapse = TRUE,
    n_cores = parallel::detectCores() - 1) {

    if (n_cores == 0)
        n_cores <- 1

    if (comp > min(rgcca$ncomp))
        stop("Selected dimension was not associated to every blocks",
            exit_code = 113)

    cat("Binding in progress...")

    mean <- weight <- sd <- occ <- list()

    if (collapse)
        J <- seq(length(rgcca$a))
    else
        J <- i_block

    for (i in J) {

        w_bind <- parallel::mclapply(w,
            function(x)
                x[[i]][, comp],
            mc.cores = n_cores)

        weight[[i]] <- rgcca$a[[i]][, comp]
        w_select <- matrix(
            unlist(w_bind),
            nrow = length(w_bind),
            ncol = length(w_bind[[1]]),
            byrow = TRUE
        )
        colnames(w_select) <- names(weight[[i]])
        rm(w_bind); gc()

        n <- seq(NCOL(w_select))

        if (is(rgcca, "sgcca")) {

            occ[[i]] <- unlist(
                parallel::mclapply(n,
                function(x)
                    sum(w_select[, x] != 0) / length(w_select[, x]),
                mc.cores = n_cores))

        }

        mean[[i]] <- unlist(parallel::mclapply(n,
            function(x) mean(w_select[,x]),
            mc.cores = n_cores
        ))
        sd[[i]] <- unlist(
            parallel::mclapply(n,
                function(x) sd(w_select[,x]),
                mc.cores = n_cores
        ))

        rm(w_select); gc()
    }

    rm(w); gc()

    occ <- unlist(occ)
    mean <- unlist(mean)
    weight <- unlist(weight)
    sd <- unlist(sd)

    cat("OK.\n", append = TRUE)

    p.vals <- pnorm(0, mean = abs(mean), sd = sd)
    tail <- qnorm(1 - .05 / 2)

    df <- data.frame(
        mean = mean,
        rgcca = weight,
        intneg = mean - tail * sd,
        intpos = mean + tail * sd,
        br = abs(mean) / sd,
        p.vals,
        BH = p.adjust(p.vals, method = "BH")
    )

    if (is(rgcca, "sgcca")) {
        index <- 8
        df$occ <- occ
    }else{
        index <- 5
        df$sign <- rep("", NROW(df))

        for (i in seq(NROW(df)))
            if (df$intneg[i]/df$intpos[i] > 0)
                df$sign[i] <- "*"

    }

    if (collapse)
        df$color <- as.factor(get_bloc_var(rgcca$a, collapse = collapse))

    zero_var <- which(df[, 1] == 0)
    if (length(zero_var) != 0)
        df <- df[-zero_var, ]

    b <- data.frame(order_df(df, index, allCol = TRUE), order = NROW(df):1)
    attributes(b)$indexes <-
        list(
            mean = "Mean bootstrap weights",
            br = "Bootstrap-ratio",
            sign = "Significant 95% interval",
            occ = "Non-zero occurences"
        )
    attributes(b)$type <- class(rgcca)
    class(b) <- c(class(b), "bootstrap")

    return(b)
}
