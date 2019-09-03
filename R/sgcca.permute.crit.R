#' Run through a set of constraint parameters c1s to select the best with permutation
#' Only one component per block for the time being
#' 
#' 
#' @param A A list that contains the J blocks of variables \mathbf{X_1}, \mathbf{X_2}, ..., \mathbf{X_J}
#' @param c1s A matrix containing sets of constraint variables, one row by set. If null, sgcca.permute takes 10 sets between min values ($1/sqrt(ncol)$) and 1
#' @param nperm Number of permutation tested for each set of constraint
#' @param C A design matrix that describes the relationships between blocks (default: complete design)
#' @param scheme The value is "horst", "factorial" or "centroid" (default: "centroid")
#' @param plot A logical, should a plot of coeffi
#' @param n_cores For linux and MacOS number of cores used for parallelisation
#' @param ncomp Number of component computed for each block
#' @value A list containing :
#' @value \item {pval}
#' @value \item {zstat}
#' @value \item {bestpenalties}
#' @value \item {permcrit} 
#' @value \item {crit}
#' @examples
#' library(RGCCA)
#' library(PMA)
#' 
#' # Create correlated blocks
#' 
#' u <- matrix(rnorm(50),ncol=1)
#' 
#' v1 <- matrix(c(rep(.5,25),rep(0,75)),ncol=1)
#' v2 <- matrix(c(rep(1,25),rep(0,25)),ncol=1)
#' v3 <- matrix(c(rep(.5,25),rep(0,175)),ncol=1)
#'
#' x1 <- u%*%t(v1) + matrix(rnorm(50*100),ncol=100)
#' x2 <- u%*%t(v2) + matrix(rnorm(50*50),ncol=50)
#' x3 <- u%*%t(v3) + matrix(rnorm(50*200),ncol=200)
#'
#' xlist <- list(x1, x2, x3)
#' 
#' # Search for the best sparsity penalities for SGCCA
#' 
#' perm.sgcca = sgcca.permute.crit(xlist, nperm = 10)
#' out.sgcca = sgcca(xlist, c1 = perm.sgcca$bestpenalties)
#' 
#' par(mfrow=c(3,1), mar = rep(2,4))
#' PlotCGH(out.sgcca$a[[1]], chrom=rep(1,ncol(x1)))
#' PlotCGH(out.sgcca$a[[2]], chrom=rep(2,ncol(x2)))
#' PlotCGH(out.sgcca$a[[3]], chrom=rep(3,ncol(x3)))
#' title(main = "SGCCA.permute, C by default", outer = TRUE, line = -2)
#' 
#'@export sgcca.permute.crit

sgcca.permute.crit <- function(
    A,
    c1s = NULL,
    nperm = 20,
    C = 1 - diag(length(A)),
    ncomp = rep(1, length(A)),
    scheme = "factorial",
    plot = FALSE,
    tol = .Machine$double.eps,
    n_cores = parallel::detectCores() - 1,
    scale = TRUE) {

    if (is.null(c1s)) {
        c1s <- matrix(NA, nrow = 10, ncol = length(A))
        for (k in 1:length(A))
            c1s[, k] <- pmin(1 / (seq(0.9, 0.01, len = 10) * sqrt(ncol(A[[k]]))), 1)
    }

    crit <- crits <- rep(0,NROW(c1s))

    for (i in 1:NROW(c1s)) {
        write(i, file = "permut.log", append = TRUE)
        out <- sgcca(
                A,
                C,
                c1 = c1s[i,],
                ncomp = ncomp,
                scheme = scheme,
                tol = tol,
                scale = scale
            )
        crits[i] <- mean(unlist(lapply(out$crit, function(x) x[length(x)])))
    }

    if (Sys.info()["sysname"] == "Windows") {

        n_cores <- parallel::detectCores() - 1
        e <- environment()
        cl <- parallel::makeCluster(n_cores)

        parallel::clusterExport(
            cl,
            c(
                "A",
                "c1s",
                "nperm",
                "C",
                "ncomp",
                "scheme",
                "out",
                "crit",
                "crits",
                "tol"
            ),
            envir = e
        )

        parallel::clusterEvalQ(cl, library(devtools))
        
        tryCatch({
            parallel::clusterEvalQ(cl, load_all("RGCCA/R/."))
        }, error = function(e) {
            warning("error : probably an issue with the localisation of RGCCA functions")
        })

        # Close cluster even if there is an error or a warning with sgcca.crit
        permcrit <- tryCatch({
            parallel::parSapply(cl, 1:nperm, function(x)
                sgcca.crit(
                    A = A,
                    C = C,
                    c1s = c1s,
                    ncomp = ncomp,
                    scheme = scheme,
                    tol = tol,
                    crit = crit,
                    scale = scale
                ))
        }, error = function(e) {
            warning("an error occured with sgcca.crit")
            return(NULL)
        })

        parallel::stopCluster(cl)
        
        if (is.null(permcrit))
            return(NULL)

    } else {
        permcrit <- simplify2array(parallel::mclapply(1:nperm,
            function(x){
                res <- sgcca.crit(
                    A = A,
                    C = C,
                    c1s = c1s,
                    ncomp = ncomp,
                    scheme = scheme,
                    tol = tol,
                    crit = crit,
                    scale = scale
                )
                return(res)
                },
            mc.cores = n_cores))
    }
                    
    pvals <- zs <- matrix(NA, nrow = NROW(c1s), ncol = NCOL(c1s) + 1)

    # chaque ligne : [critère 1, critère 2, valeur pvals ou zs]
    for (i in 1:NROW(c1s)) {
        pvals[i,] <- c(c1s[i,], mean(permcrit[i,] >= crits[i]))
        zs[i,] <- c(c1s[i,], (crits[i] - mean(permcrit[i,])) / (sd(permcrit[i,])))  # + 0.05))
    }

    bestpenalties <- c1s[which.max(zs[, length(A) + 1]), 1:length(A)]

    sgcca.best <- sgcca(
            A,
            C = C,
            c1 = bestpenalties,
            ncomp = ncomp,
            scheme = scheme
        )

    out <- list(
            pvals = pvals,
            zstat = zs,
            bestpenalties = bestpenalties,
            sgcca.best = sgcca.best,
            permcrit = permcrit,
            crit = crits,
            penalties = c1s
        )

    if (plot) {

        par(mfrow = c(2, 1), mar = c(2, 4, 2, 1))

        plot(1:NROW(c1s),
            out$crit,
            ylim = c(0, max(crits, permcrit)),
            ylab = "pvals")

        for (i in 1:nperm)
            points(1:NROW(c1s), out$permcrit[, i], col = "green")

        for (j in 1:(NROW(c1s) - 1))
            segments(j, out$crit[j], j + 1, out$crit[j + 1])

        plot(1:NROW(c1s), out$zstat[, length(A) + 1], ylab = "zstat")
        for (j in 1:(NROW(c1s) - 1)) 
            segments(
                    j, 
                    out$zstat[j, length(A) + 1],
                    j + 1, 
                    out$zstat[j + 1, length(A) + 1]
                )

    }

    return(out)
}
