#' Extract statistics from a fitted bootstrap object
#'
#' This function extracts statistical information from a fitted bootstrap
#' object (see \code{\link[RGCCA]{bootstrap}}).
#'
#' @inheritParams bootstrap
#' @inheritParams plot_histogram
#' @inheritParams plot_var_2D
#' @inheritParams plot.rgcca
#' @inheritParams plot_var_1D
#' @param bars Character string giving the variability among "sd" (error bar of
#' standard deviation) , "stderr" (error bar of standard deviation divided by 
#' sqrt(n)) or "quantile" (error bar of 0.05/0.95-quantiles)
#' @param b A fitted bootstrap object (see  \code{\link[RGCCA]{bootstrap}} )
#' @param display_order A logical value to display the order of the variables
#' @param adj.method character string indicating the method used for p-value 
#' adjustment (default: fdr - a.k.a Benjamini-Hochberg correction)
#' @return A dataframe containing:
#' \itemize{
#' \item 'mean' for the mean of the bootstrap weights (non-null for SGCCA)
#' \item 'estimate' for RGCCA weights
#' \item 'sd' for the standard error of the (non-null in case of SGCCA) 
#' bootstrap weights
#' \item 'lower/upper_bound' for the lower and upper intervals calculated 
#' according to the 'bars' parameter
#' \item 'bootstrap_ratio' for the mean of the bootstrap weights / their 
#' standard error
#' \item 'pval' for p-value (see details)
#' \item 'adjust.pval' for adjusted p-value (default value: fdr (Benjamini-
#' Hochberg correction))
#' \item 'occurrences' for non-zero occurences (for SGCCA) 
#' \item 'sign' for significant ('*') or not ('ns') p-value (alpha = 0.05) 
#' (see details)
#' }
#' @details 
#' For RGCCA, the p-value is computed by assuming that the ratio of the blocks 
#' weight values to the bootstrap estimate of the standard deviation follows 
#' the standardized normal distribution. 
#' By including sparsity (with "sgcca","spls" or "spca"), the frequency of a 
#' selected variable may depend on both the level of sparsity and the total 
#' number of variables in each block. 
#' For a random selection of the variable within the block, the number of 
#' occurrences (0 or 1) follows a Bernoulli distribution with the parameter 
#' p = proportion of selected variables in the block. 
#' This proportion is estimated by the average number of selected variables 
#' over all bootstraps divided by the total number of variables in each block 
#' (p_j). On a larger number of bootstrap samples, the number of occurrences 
#' follows a binomial distribution B(n,p) with n=number of bootstraps. 
#' The test is based on the following null hypothesis: "the variable is 
#' randomly selected according to B(n,p)". This hypothesis is rejected when 
#' the number of occurrences is higher than the 1-(0.05/p_j)th quantile 
#' @examples
#' # Bootstrap confidence intervals and p-values for RGCCA 
#' data("Russett")
#' blocks = list(agriculture = Russett[, seq(3)], 
#'               industry = Russett[, 4:5], 
#'               politic = Russett[, 6:11])
#'               
#' rgcca_out = rgcca(blocks)
#' boot = bootstrap(rgcca_out, 5, n_cores = 1)
#' get_bootstrap(boot)
#' 
#' # Stability of the selected variables for SGCCA
#'  
#' @export
#' @importFrom stats pt pbinom
#' @seealso \code{\link[RGCCA]{bootstrap}}, 
#' \code{\link[RGCCA]{plot.bootstrap}}, 
#' \code{\link[RGCCA]{print.bootstrap}} 
get_bootstrap <- function(
    b,
    comp = 1,
    block = length(b$bootstrap[[1]]),
    bars="quantile",
    display_order=TRUE,
    adj.method = "fdr") {

    stopifnot(is(b, "bootstrap"))
    check_ncol(b$rgcca$Y, block)
    check_blockx("block", block, b$rgcca$call$blocks)
    check_compx("comp", comp, b$rgcca$call$ncomp, block)
    match.arg(bars,c("quantile", "sd", "stderr"))

    bootstrapped=b$bootstrap[[comp]][[block]]
    n_boot=sum(!is.na(bootstrapped[1, ]))
    if(tolower(b$rgcca$call$type) %in% c("spls", "spca", "sgcca"))
    {
      mean = apply(bootstrapped, 1, 
                   function(x)
                     ifelse(sum(x!=0)>=1, mean(x[x != 0], na.rm = T), 0))
        
      sd = apply(bootstrapped, 1, 
                   function(x)
                     ifelse(sum(x!=0)>1, sd(x[x != 0], na.rm = T), 0))
    }
    else
    {
        mean = apply(bootstrapped, 1, function(x) mean(x, na.rm = T))
        sd = apply(bootstrapped, 1, function(x) sd(x, na.rm = T))
    }
    weight <- b$rgcca$a[[block]][, comp]
    occ <- apply(bootstrapped, 1,
            function(x)
                sum(x!= 0, na.rm=T) )
    
     # p values if occurrences
    if(tolower(b$rgcca$call$type) %in% c("spls", "spca", "sgcca"))
    {
        n_boot=ifelse(!is.null(dim(b[[1]][[1]][[1]])), 
                      dim(b[[1]][[1]][[1]])[2],
                      length(b[[1]][[1]][[1]]))
        nvar = length(b$bootstrap[[1]][[block]][, 1])
        avg_n_occ = sum(occ)/n_boot
        probComp = avg_n_occ/nvar
        q1 = qbinom(size=n_boot, prob = probComp, p=1-0.05/nvar)
        p.vals = pbinom(occ, size = n_boot, prob = probComp, lower.tail = FALSE)
    }
    else
    {
        tail <- qnorm(1 - .05 / 2)
        p.vals <- 2 * pnorm(abs(weight)/sd, lower.tail = FALSE)
    }
  
   
    if(bars=="quantile")
    {
        lower_bound=apply(bootstrapped, 1, 
                         function(y){return(quantile(y, 0.025))})
        upper_bound=apply(bootstrapped, 1,
                         function(y){return(quantile(y, 0.975))})
    }
    if(bars=="sd")
    {
        length_bar = sd
        lower_bound = mean-length_bar
        upper_bound = mean+length_bar
     }
    if(bars=="stderr")
    {
        length_bar = sd/sqrt(n_boot)
        lower_bound = mean-length_bar
        upper_bound = mean+length_bar
    }
    
    df <- data.frame(
        mean = mean,
        estimate = weight,
        sd = sd,
        lower_bound = lower_bound,
        upper_bound = upper_bound,
        bootstrap_ratio = abs(mean) / sd,
        pval = p.vals,
        adjust.pval = p.adjust(p.vals, method = adj.method)
    )

    if (tolower(b$rgcca$call$type) %in% c("spls", "spca", "sgcca")){
      df$occurrences <- occ
      index <- which(colnames(df)=="occurrences")
      db <- data.frame(order_df(df, index, allCol = TRUE), order = NROW(df):1) 
      if(!display_order){
        db <- data.frame(order_df(df, index, allCol = TRUE)   )
        db <- db[, c("occurrences", "mean", "estimate", "sd",
                     "lower_bound", "upper_bound", 
                     "pval", "adjust.pval")] 
      }
      if(display_order)
            db <- data.frame(order_df(df, index, allCol = TRUE), 
                             order = NROW(df):1)
    }
    else{
      df$sign <- rep("NS", NROW(df))
      for (i in seq(NROW(df)))
            if(p.vals[i]<0.05) df$sign[i] <- "*"
      if(!display_order){
        index <- which(colnames(df)=="mean")
        db <- data.frame(order_df(df, index, allCol = TRUE)   )
        db <- db[, c("mean", "estimate", "sd", 
                     "lower_bound", "upper_bound",
                     "pval", "adjust.pval")] 
      }
      if(display_order){
        index <- which(colnames(df) == "mean")
        db <- data.frame(order_df(df, index, allCol = TRUE), 
                             order = NROW(df):1) 
      }
    }
    
    zero_var <- which(df[, 1] == 0)
    if (NROW(df) > 1 && length(zero_var) != 0) db <- db[-zero_var, ]
   
    attributes(db)$indexes <- list(
            estimate = "RGCCA weights",
            bootstrap_ratio = "Bootstrap-ratio",
            sign = "Significance",
            occurrences = "Non-zero occurrences",
            mean = "Mean bootstrap weights")
    
    attributes(db)$type <- class(b$rgcca)
    attributes(db)$n_boot <- n_boot
    attributes(db)$n_blocks <- length(b$rgcca$call$blocks)
    class(db) <- c(class(db), "df_bootstrap")
    return(db)
}
