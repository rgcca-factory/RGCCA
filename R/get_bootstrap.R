#' Extract a bootstrap
#'
#' Extract statistical information from a bootstrap
#'
#' @inheritParams bootstrap
#' @inheritParams plot_histogram
#' @inheritParams plot_var_2D
#' @inheritParams plot_var_1D
#' @param bars A character among "sd" for standard deviations, "stderr" for the standard error, "ci" for confidence interval of scores and "cim" for the confidence intervall of the mean.
#' @param b A list of list weights (one per bootstrap per blocks)
#' @param display_order If TRUE the order is indicated
#' @return A matrix containing the means, 95% intervals, bootstrap ratio and p-values
#' @examples
#' library(RGCCA)
#' data("Russett")
#' blocks = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
#'     politic = Russett[, 6:11] )
#' rgcca_out = rgcca(blocks)
#' boot = bootstrap(rgcca_out, 5, superblock = FALSE, n_cores = 1)
#' get_bootstrap(boot, n_cores = 1)
#' @export
#' @importFrom stats pt
get_bootstrap <- function(
    b,
    comp = 1,
    i_block = length(b$bootstrap[[1]]),
    bars="sd",
    collapse = FALSE,
    n_cores = parallel::detectCores() - 1,
    display_order=FALSE)
    {
    stopifnot(is(b, "bootstrap"))

    check_compx("comp", comp, b$rgcca$call$ncomp, i_block)
    check_ncol(b$rgcca$Y, i_block)
    check_blockx("i_block", i_block, b$rgcca$call$blocks)
    check_boolean("collapse", collapse)
    check_integer("n_cores", n_cores, 0)
    bootstrapped=b$bootstrap[[comp]][[i_block]]
    n_boot=dim(bootstrapped)[2]
    mean=apply(bootstrapped,1,mean)
    sd=apply(bootstrapped,1,sd)
    weight <- b$rgcca$a[[i_block]][, comp]
    occ <- apply(bootstrapped,1,
            function(x)
                sum(x!= 0) / length(x))
    print("in")
    p.vals <- 2 * pt(abs(weight)/sd, lower.tail = FALSE, df = n_boot - 1)
    tail <- qt(1 - .05 / 2, df = n_boot - 1)
    
    if(bars=="sd")
    {
        length_bar=sd
     }
    if(bars=="stderr")
    {
        length_bar=sd/sqrt(n_boot)
    }
    if(bars=="ci")
    {
        length_bar=tail*sd
    }
    # if(bars=="cim")
    # {
    #     length_bar=tail*sd/sqrt(n_boot)
    # }
    # 
    df <- data.frame(
        mean = mean,
        estimate = weight,
        sd=sd,
        lower_band = mean -  length_bar,
        upper_band = mean +  length_bar,
        bootstrap_ratio = abs(mean) / sd,
        p.vals,
        BH = p.adjust(p.vals, method = "BH")
    )
print("before if")
    if (tolower(b$rgcca$call$type) %in% c("spls", "spca", "sgcca")) {
          df$occurrences <- occ
          index <- which(colnames(df)=="occurrences")
          print("here")
        db <- data.frame(order_df(df, index, allCol = TRUE), order = NROW(df):1) 
        print("out")
    }
    else
    {
        df$sign <- rep("", NROW(df))
        
        for (i in seq(NROW(df)))
            #if (df$lower_band[i]/df$upper_band[i] > 0)
            #if (abs(df$mean[i]/(df$sd[i]/sqrt(n_boot))) > qt(0.95/2,df=n_boot-1))
            if(p.vals[i]<0.05)
                df$sign[i] <- "*"
        if(!display_order)
        {
            df=df[,c("mean","estimate","sd","lower_band","upper_band","p.vals","BH")] 
            index <- which(colnames(df)=="mean")
            db <- data.frame(order_df(df, index, allCol = TRUE)   )
        }
        if(display_order)
        {
            index <- which(colnames(df)=="mean")
            db <- data.frame(order_df(df, index, allCol = TRUE), order = NROW(df):1) 
        }     
    }
    
   # rm(b); gc()

    if (collapse)
        db$color <- as.factor(get_bloc_var(b$rgcca$a, collapse = collapse))

    zero_var <- which(df[, 1] == 0)
    if (NROW(df) > 1 && length(zero_var) != 0)
        db <- db[-zero_var, ]
   
       attributes(db)$indexes <-
        list(
            estimate = "RGCCA weights",
            bootstrap_ratio = "Bootstrap-ratio",
            sign = "Significant 95% \nbootstrap interval",
            occurrences = "Non-zero occurences",
            mean = "Mean bootstrap weights"
        )
    attributes(db)$type <- class(rgcca)
    attributes(db)$n_boot <- n_boot
    class(db) <- c(class(db), "df_bootstrap")
    return(db)
}
