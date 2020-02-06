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
#' rgcca_out = rgcca(blocks)
#' boot = bootstrap(rgcca_out, 2, FALSE, n_cores = 1)
#' get_bootstrap(rgcca_out, boot, n_cores = 1)
#' @export
#' @importFrom stats pt
get_bootstrap <- function(
   
    w,
    comp = 1,
    i_block = length(w[[1]]),
    collapse = FALSE,
    n_cores = parallel::detectCores() - 1) {

    check_ncol(w$w$rgcca$Y, i_block)

    if (n_cores == 0)
        n_cores <- 1

    if (collapse && w$rgcca$call$superblock) {
        w$rgcca$a <- w$rgcca$a[-length(w$rgcca$a)]
        if (i_block > length(w$rgcca$a))
            i_block <- length(w$rgcca$a)
    }

    if (comp > min(w$rgcca$call$ncomp))
        stop("Selected dimension was not associated to every blocks",
             exit_code = 113)

    cat("Binding in progress...")

    mean <- weight <- sd <- occ <- list()

    if (collapse)
        J <- seq(length(w$rgcca$a))
    else
        J <- i_block

    for (i in J) {

        w_bind <- parallel::mclapply(
            w,
            function(x) x[[i]][, comp],
             mc.cores = n_cores)

        weight[[i]] <- w$rgcca$a[[i]][, comp]
        w_select <- matrix(
            unlist(w_bind),
            nrow = length(w_bind),
            ncol = length(w_bind[[1]]),
            byrow = TRUE
        )
        colnames(w_select) <- names(weight[[i]])
        rm(w_bind); gc()

        n <- seq(NCOL(w_select))

        if (w$rgcca$call$type %in% c("spls", "spca", "sgcca")) {

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
        if (NCOL(w_select) == 1)
            sd[[i]] <- 1
        else
            sd[[i]] <- unlist(
                parallel::mclapply(n,
                                   function(x) sd(w_select[,x]),
                                   mc.cores = n_cores
                ))

        rm(w_select); gc()
    }
    
    n_boot <- length(w)
    rm(w); gc()

    occ <- unlist(occ)
    mean <- unlist(mean)
    weight <- unlist(weight)
    sd <- unlist(sd) / sqrt(n_boot)

    cat("OK.\n", append = TRUE)

    p.vals <- 2 * pt(abs(weight)/sd, lower.tail = FALSE, df = n_boot - 1)
    tail <- qt(1 - .05 / 2, df = n_boot - 1)

#

    df <- data.frame(
        mean = mean,
        estimate = weight,
        lower_band = mean - (tail * sd/sqrt(n_boot)),
        upper_band = mean + (tail * sd/sqrt(n_boot)),
        bootstrap_ratio = abs(mean) / sd,
        p.vals,
        BH = p.adjust(p.vals, method = "BH")
    )

    if (w$rgcca$call$type %in% c("spls", "spca", "sgcca")) {
        index <- 8
        df$occurrences <- occ
    }else{
        index <- 5
        df$sign <- rep("", NROW(df))
        
        for (i in seq(NROW(df)))
            if (df$intneg[i]/df$intpos[i] > 0)
                df$sign[i] <- "*"
        
    }

    if (collapse)
        df$color <- as.factor(get_bloc_var(w$rgcca$a, collapse = collapse))

    zero_var <- which(df[, 1] == 0)
    if (NROW(df) > 1 && length(zero_var) != 0)
        df <- df[-zero_var, ]

    b <- data.frame(order_df(df, index, allCol = TRUE), order = NROW(df):1)
    attributes(b)$indexes <-
        list(
            mean = "Mean bootstrap weights",
            bootstrap_ratio = "Bootstrap-ratio",
            sign = "Significant 95% \nbootstrap interval",
            occurrences = "Non-zero occurrences",
            estimate = "RGCCA weights"
        )
    attributes(b)$type <- class(rgcca)
    attributes(b)$n_boot <- n_boot
    class(b) <- c(class(b), "bootstrap")

    return(b)
}
