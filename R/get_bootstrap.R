#' Extract a bootstrap
#'
#' Extract statistical information from a bootstrap
#'
#' @inheritParams bootstrap
#' @inheritParams plot_histogram
#' @inheritParams plot_var_2D
#' @inheritParams plot_var_1D
#' @param bars A character among "sd" for standard deviations, "stderr" for the standard error, "ci" for confidence interval of scores and "cim" for the confidence intervall of the mean.
#' @param b A bootstrap object (see  \code{\link[RGCCA]{bootstrap}} )
#' @param display_order A logical value to display the order of the variables
#' @return A matrix containing the means, 95\% intervals, bootstrap ratio, p-values and other statistics (see details)
#' @details 
#' \itemize{
#' \item 'mean' for the mean of the bootstrap weights
#' \item 'estimate' for RGCCA weights
#' \item 'sd' for the standard error of the bootstrap weights
#' \item 'lower/upper_band' for the lower and upper intervals from to the 'bar' parameter
#' \item 'bootstrap_ratio' for the mean of the bootstrap weights / their standard error
#' \item 'p.vals' for p-values. In the case of SGCCA, the occurrences of the 
#' bootstrap weights are distributed according to the binomial distribution 
#' while for RGCCA, their absolute value weighted by their standard deviation 
#' follows a normal distribution.
#' \item 'BH' for Benjamini-Hochberg p-value adjustments
#' \item 'occurrences' for non-zero occurences (for SGCCA) 
#' \item 'sign' for significant ('*') or not ('ns') p-values (alpha = 0.05)
#' }
#' @examples
#' library(RGCCA)
#' data("Russett")
#' blocks = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
#'     politic = Russett[, 6:11] )
#' rgcca_out = rgcca(blocks)
#' boot = bootstrap(rgcca_out, 5, n_cores = 1)
#' get_bootstrap(boot, n_cores = 1)
#' @export
#' @importFrom stats pt pbinom
#' @seealso \code{\link[RGCCA]{bootstrap}}, \code{\link[RGCCA]{plot.bootstrap}} , \code{\link[RGCCA]{print.bootstrap}} 
get_bootstrap <- function(
    b,
    comp = 1,
    i_block = length(b$bootstrap[[1]]),
    bars="sd",
    collapse = FALSE,
    n_cores = parallel::detectCores() - 1,
    display_order=TRUE)
    {
    stopifnot(is(b, "bootstrap"))

    check_compx("comp", comp, b$rgcca$call$ncomp, i_block)
    check_ncol(b$rgcca$Y, i_block)
    check_blockx("i_block", i_block, b$rgcca$call$blocks)
    check_boolean("collapse", collapse)
    check_integer("n_cores", n_cores, 0)
    bootstrapped=b$bootstrap[[comp]][[i_block]]
    #n_boot=dim(bootstrapped)[2]
    n_boot=sum(!is.na(bootstrapped[1,]))
    if(tolower(b$rgcca$call$type) %in% c("spls", "spca", "sgcca"))
    {
        mean=apply(bootstrapped,1,function(x){mean(x[x!=0],na.rm=T)})
        sd=apply(bootstrapped,1,function(x){sd(x[x!=0],na.rm=T)})
    }
    else
    {
        mean=apply(bootstrapped,1,function(x){mean(x,na.rm=T)})
        sd=apply(bootstrapped,1,function(x){sd(x,na.rm=T)})
    }
    weight <- b$rgcca$a[[i_block]][, comp]
    occ <- apply(bootstrapped,1,
            function(x)
                sum(x!= 0,na.rm=T) )
    
     # p values if occurrences
    if(tolower(b$rgcca$call$type) %in% c("spls", "spca", "sgcca"))
    {
        n_boot=ifelse(!is.null(dim(b[[1]][[1]][[1]])),dim(b[[1]][[1]][[1]])[2],length(b[[1]][[1]][[1]]))
        nvar=length(b$bootstrap[[1]][[i_block]][,1])
        avg_n_occ=sum(occ)/n_boot
        probComp= avg_n_occ/nvar
        q1=qbinom(size=n_boot,prob=probComp,p=1-0.05/nvar)
        p.vals=pbinom(occ,size=n_boot,prob=probComp,lower.tail=FALSE)
        
    }
    else
    {
        tail <- qt(1 - .05 / 2, df = n_boot - 1)
        p.vals <- 2 * pt(abs(weight)/sd, lower.tail = FALSE, df = n_boot - 1)
    }
  
   
    if(bars=="quantile")
    {
        lower_band=apply(bootstrapped,1,function(y){return(quantile(y,0.05))})
        upper_band=apply(bootstrapped,1,function(y){return(quantile(y,0.95))})
    }
    if(bars=="sd")
    {
        length_bar=sd
        lower_band=mean-length_bar
        upper_band=mean+length_bar
     }
    if(bars=="stderr")
    {
        length_bar=sd/sqrt(n_boot)
        lower_band=mean-length_bar
        upper_band=mean+length_bar
    }
    if(bars=="ci")
    {
        length_bar=tail*sd
        lower_band=mean-length_bar
        upper_band=mean+length_bar
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
        lower_band = lower_band,
        upper_band = upper_band,
        bootstrap_ratio = abs(mean) / sd,
        p.vals,
        BH = p.adjust(p.vals, method = "BH")
    )

    if (tolower(b$rgcca$call$type) %in% c("spls", "spca", "sgcca")) 
    {
          df$occurrences <- occ
          index <- which(colnames(df)=="occurrences")
        db <- data.frame(order_df(df, index, allCol = TRUE), order = NROW(df):1) 
        if(!display_order)
        {
             db <- data.frame(order_df(df, index, allCol = TRUE)   )
             db <- db[,c("occurrences","mean","estimate","sd","lower_band","upper_band","p.vals","BH")] 
        }
        if(display_order)
        {
            db <- data.frame(order_df(df, index, allCol = TRUE), order = NROW(df):1) 
        }  
   
    }
    else
    {
        df$sign <- rep("NS", NROW(df))
        for (i in seq(NROW(df)))
            #if (df$lower_band[i]/df$upper_band[i] > 0)
            #if (abs(df$mean[i]/(df$sd[i]/sqrt(n_boot))) > qt(0.95/2,df=n_boot-1))
            if(p.vals[i]<0.05)
                df$sign[i] <- "*"
        if(!display_order)
        {
             index <- which(colnames(df)=="mean")
            db <- data.frame(order_df(df, index, allCol = TRUE)   )
            db <- db[,c("mean","estimate","sd","lower_band","upper_band","p.vals","BH")] 
            
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
    attributes(db)$type <- class(b$rgcca)
    attributes(db)$n_boot <- n_boot
    attributes(db)$n_blocks <- length(b$rgcca$call$blocks)
    class(db) <- c(class(db), "df_bootstrap")
    return(db)
}
