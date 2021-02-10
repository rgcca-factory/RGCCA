# Other checking functions
#-----------------------------
check_blockx <- function(x, y, blocks){
    x <- check_min_integer(x, y, blocks)
    if (y > length(blocks))
        stop_rgcca(
            paste0(
                x,
                " should be lower than ",
                length(blocks),
                " (that is the number of blocks)."
            ),
            exit_code = 133
        )
    return(x)
}

check_boolean <- function(x, y = x, type = "scalar") {
    
    if (is.null(y))
        y <- x
    
    if (type == "scalar")
        x = ""
    
    if (any(is.na(y)))
        stop_rgcca(paste(x, "should not be NA."))
    
    if (!is(y, "logical"))
        stop_rgcca(paste(x, "should be TRUE or FALSE"))
    
    if (type == "scalar" && length(y) != 1)
        stop_rgcca(paste(x, "should be of length 1."))
}

check_colors <- function(colors, n = 3){
    
    if (!is.null(colors)){
        colors <- as.vector(colors)
        for (i in colors) {
            if (!is.na(i) && !(i  %in%  colors()) && is.character2(i)
                && regexpr("^#{1}[a-zA-Z0-9]{6,8}$", i) < 1)
                stop_rgcca("colors must be in colors() or a rgb character.")
        }
    }
}

check_compx <- function(x, y, ncomp, blockx) {
    res <- check_integer(x, y, min = 1)
    if (y > ncomp[blockx]) {
        stop_rgcca(
            paste0(
                x,
                " is equals to ",
                y,
                " and should be less than or equal to ",
                ncomp[blockx],
                " (the number of block component)."
            ),
            exit_code = 128
        )
    }
 return(res)
}

# Check the format of the connection matrix
#
# @inheritParams rgccad
# @inheritParams set_connection
check_connection <- function(C, blocks) {
    
    msg <- "The design matrix C should"
    
    if (!isSymmetric.matrix(unname(C)))
        stop_rgcca(paste(msg, "be symmetric."), exit_code = 103)
    
    # d <- unique(diag(C))
    # if (length(d) != 1 || d != 0)
    #     stop_rgcca("The diagonal of the connection matrix file should be 0.",
    #         exit_code = 105)
    
    x <- C >=0 & C<=1
    if (sum(!x)!=0)
        stop_rgcca(paste(msg, "contain numbers between 0 or 1."), exit_code = 106)
    
    if (all(C == 0))
        stop_rgcca(paste(msg, "not contain only 0."), exit_code = 107)
    
    invisible(check_size_blocks(blocks, "connection matrix", C))
    
    if(is.null(rownames(C)) || is.null(colnames(C)))
        rownames(C) <- names(blocks) -> colnames(C)
    
    if (!all(rownames(C) %in% names(blocks)) || 
        !all(colnames(C) %in% names(blocks)))
        stop_rgcca(paste(msg,
                         "have the rownames and the colnames that match with 
                         the names of the blocks."),
                   exit_code = 108)
    
    return(C)
    # TODO: warning if superblock = TRUE
}
check_file <- function(f) {
    # Check the existence of a path f: A character giving the path of a file
    
    if (!file.exists(f))
        stop_rgcca(paste(f, "does not exist."), exit_code = 101)
    
}
# y <- runif(6)
# check_integer("y", y, "vector", T, min = 0)
# y <- matrix(runif(6, 1, 10), 2, 3)
# check_integer("y", y, "matrix", T)
# check_integer("y", y, "matrix")
# check_integer(NA)
# check_integer(.1)
# check_integer(c(1:2))
# check_integer("x", c(0, 0), type = "vector")
check_integer <- function(x, y = x, type = "scalar", float = FALSE, min = 1) {
    
    if (is.null(y))
        y <- x
    
    if (type %in% c("matrix", "data.frame"))
        y_temp <- y
    
    y <- tryCatch(
        as.double(as.matrix(y)),
        warning = function(w)
            stop_rgcca(paste(x, "should be numeric."))
    )
    
    if (any(is.na(y)))
        stop_rgcca(paste(x, "should not be NA."))
    
    if (!is(y, "numeric"))
        stop_rgcca(paste(x, "should be numeric."))
    
    if (type == "scalar" && length(y) != 1)
        stop_rgcca(paste(x, "should be of length 1."))
    
    if (!float)
        if (any((y %% 1) != 0)) {
            stop_rgcca(paste(x, "should be an integer."))
        }
        y = as.integer(y)
    
    if (all(y < min)) # TODO: Ask why it is all and not any
        stop_rgcca(paste0(x, " should be higher than or equal to ", min, "."))
    
    if (type %in% c("matrix", "data.frame"))
        y <- matrix(
            y, 
            dim(y_temp)[1], 
            dim(y_temp)[2],
            dimnames = dimnames(y_temp)
        )
    
    if (type == "data.frame")
        as.data.frame(y)
    
    return(y)
}
check_lower_blocks <- function(x, y, blocks)
    if (y > length(blocks))
        stop_rgcca(
            paste0(
                x,
                " should be lower than ",
                length(blocks),
                " (the maximum number of blocks), not be equal to ",
                y,
                "."
            ),
            exit_code = 133
        )
check_method <- function(method) {
    analysis <- c("rgcca", "sgcca", "mgcca", "pca", "spca", "pls", "spls", 
      "cca", "ifa", "ra", "gcca", "maxvar", "maxvar-b", 
      "maxvar-a", "mcoa","cpca-1", "cpca-2", "cpca-4", 
      "hpca", "maxbet-b", "maxbet", "maxdiff-b", "maxdiff", 
      "maxvar-a", "sabscor", "ssqcor", "ssqcor", "ssqcov-1", 
      "ssqcov-2", "ssqcov", "sumcor", "sumcov-1", "sumcov-2", 
      "sumcov", "sabscov", "sabscov-1", "sabscov-2")
      
    if (!tolower(method) %in% analysis)
        stop_rgcca(
            paste0("Please select one type among the following 
            type: ", paste(analysis, collapse = ", ")),
            exit_code = 112
        )
}

check_min_integer <- function(x , y, blocks, min = 1) {
    x <- check_integer(x, y, min = min)
    check_lower_blocks(x, y, blocks)
    return(x)
}
check_nblocks <- function(blocks, type) {
    if (tolower(type) == "pca") {
        msg <- "Only one block is"
        exit_code <- 110
    } else{
        msg <- "Two blocks are"
        exit_code <- 111
    }
    
    stop_rgcca(
        paste0(
            length(blocks),
            " blocks used in the analysis. ",
            msg ,
            " required for a ",
            type,
            "."
        ),
        exit_code = exit_code
    )
}
# If less than 2 columns, do not run
# x, a list of matrix
# i_block, the position of the tested matrix in the list
# return an error or NULL
check_ncol <- function(x, i_block) {
    if (NROW(x[[i_block]]) < 2) {
        stop_rgcca(
            "This output is available only for more than one variable."
        )
    }
}

check_ncomp <- function(ncomp, blocks, min = 1) {
    
    ncomp <- elongate_arg(ncomp, blocks)
    ncomp <- sapply(
        seq(length(ncomp)),
        function(x){
            y <- check_integer("ncomp", ncomp[x], min = min)
            if (y > NCOL(blocks[[x]])) {
                stop_rgcca(
                    paste0(
                        "ncomp[", x, "] should be comprise between ", min ," and ", NCOL(blocks[[x]]), " (that is the number of variables of block ", x, ")."
                    ),
                    exit_code = 126
                )
            }
            else
                return(y)
        }
    )
    
    check_size_blocks(blocks, "ncomp", ncomp)
    return(ncomp)
}

check_ranks <- function(ranks, blocks, min = 1) {
  ranks <- elongate_arg(ranks, blocks)
  ranks <- sapply(
    seq(length(ranks)),
    function(x){
      y <- check_integer("ranks", ranks[x], min = min)
      if (any(y > dim(blocks[[x]])[-1])) {
        stop_rgcca(
          paste0(
            "ranks[", x, "] should be comprise between ", min ,
            " and ", min(dim(blocks[[x]])), " (that is the number of
            variables of the smallest mode of block ", x, ")."
          ),
          exit_code = 126 # TODO: look up exit codes
        )
      }
      else
        return(y)
    }
  )
  
  check_size_blocks(blocks, "ranks", ranks)
  return(ranks)
}

check_reg_matrices <- function(regularisation_matrices, blocks) {
  if (is.null(regularisation_matrices)) return(regularisation_matrices)
  for (j in 1:length(blocks)) {
    DIM = dim(blocks[[j]])
    if (length(DIM) < 3) {
      message(paste0("Regularization matrices are not available for matrix 
                     blocks so regularisation_matrices[[", j, "]] has been set
                     to NULL."))
      regularisation_matrices[[j]] = NULL
    }
    if (!is.null(regularisation_matrices[[j]])) {
      reg_DIM = lapply(regularisation_matrices[[j]], dim)
      if (any(sapply(reg_DIM, function(x) x[1]) != sapply(reg_DIM, function(x) x[2]))) {
        stop_rgcca("regularisation_matrices matrices must be square matrices")
      }
      if (length(regularisation_matrices[[j]]) != length(DIM) - 1) {
        stop_rgcca(paste0("There should be as many regularisation_matrices 
                          matrices as modes in the block. Mismatch found for
                          block ", j, "."))
      }
      if (any(sapply(reg_DIM, function(x) x[1]) != DIM[-1])) {
        stop_rgcca(paste0("regularisation_matrices matrices should match the 
                          mode dimensions. Mismatch found for block ", j, "."))
      }
    }
  }
  return(regularisation_matrices)
}

# Check if a dataframe contains no qualitative variables
# @inheritParams load_blocks
# @param df A dataframe or a matrix
# @param fo A character giving the name of the tested file
# @param warn_separator A bolean to print warning for bad separator use
check_quantitative <- function(df, fo, header = FALSE, warn_separator = FALSE) {
    qualitative <- is.character2(df, warn_separator = TRUE)
    
    if (qualitative) {
        msg <- paste(
            fo,
            "contains qualitative data. Please, transform them in a disjunctive table."
        )
        
        if (!header)
            msg <- paste0(msg, "Possible mistake: header parameter is disabled, check if the file doesn't have one.")
        
        stop_rgcca(paste(msg, "\n"), exit_code = 100)
    }
    
}

check_response <- function(response = NULL, df = NULL) {
    
    if (!is.null(response)) {
        qualitative <- is.character(response)
        
        # if (length(qualitative) > 1)
        #     stop_rgcca(
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
                warning("There is multiple columns in the response block. By default, only the first column will be considered.")
                return(as.matrix(response[, 1]))
            }
        }
        
        return(response)
        
    } else
        return(rep(1, NROW(df[[1]])))
    
}

# Test on the sign of the correlation
check_sign_comp <- function(rgcca, w){
    
    w1 <- rgcca$a
    
    for (k in seq(length(w))) {
        if(NCOL(w[[k]])>1)
        {
            for (j in seq(NCOL(w[[k]]))) {
                
                res <- cor(w1[[k]][, j], w[[k]][, j])
                if (!is.na(res) && res  < 0)
                    w[[k]][, j] <- -1 * w[[k]][, j]
            }
        }
        else
        {
            res <- cor(w1[[k]], w[[k]])
            if (!is.na(res) && res  < 0)
                w[[k]] <- -1 * w[[k]]
        }
    }
    
    return(w)
}

check_size_blocks <- function(blocks, x, y = x) {
    
    if (identical(x, y))
        x <- ""
    if (any(class(y) %in% c("matrix", "data.frame"))) {
        dim_y <- NCOL(y)
        dim_type <- "number of columns"
    }else{
        dim_y <- length(y)
        dim_type <- "size"
    }
    
    if (dim_y != length(blocks))
        stop_rgcca(
            paste0(
                x,
                " should have the same ", 
                dim_type , 
                " (actually ",
                dim_y,
                ") than the number of blocks (",
                length(blocks),
                ")."
            ),
            exit_code = 130
        )
    else
        return(TRUE)
}


# Print warning if file size over
check_size_file <- function(filename) {
    size <- file.size(filename)
    if (size > 5e+06)
        # warning(paste0('The size of ', filename, ' is over 5 Mo (',
        #  round(size / 1E6, 1), ' Mo). File loading could take some times...'),
        message("File loading in progress ...")
}

check_spars <- function(blocks, tau, type = "rgcca") {
    # sparsity : A vector of integer giving the spasity parameter for SGCCA (sparsity)
    # Stop the program if at least one sparsity parameter is not in the required interval
    
    if (tolower(type) == "sgcca") {
        #the minimum value avalaible
        min_sparsity <- lapply(blocks, function(x) 1 / sqrt(NCOL(x)))
        
        # Check sparsity varying between 1/sqrt(pj) and 1
        tau <- mapply(
            function(x, y) {
                x <- check_integer("sparsity", x, float = TRUE, min = 0)
                if (x < y | x > 1)
                    stop_rgcca(
                        paste0(
                            "Sparsity parameter is equals to ",
                            x,
                            ". For SGCCA, it must be comprise between 1/sqrt(number_column) (i.e., ",
                            toString(unlist(
                                lapply(min_sparsity, function(x)
                                    ceiling(x * 100) / 100)
                            ))
                            ,
                            ") and 1."
                        ),
                        exit_code = 132
                    )
                else
                    x
            }, tau, min_sparsity)
    }
    
    invisible(tau)
}
# #' @export
check_superblock <- function(is_supervised = NULL, is_superblock = NULL, verbose = TRUE) {
    if (!is.null(is_supervised)) {
        if (verbose)
            warn_connection("supersized method with a response")
        if (is_superblock) {
            if (!is.null(is_superblock) && verbose)
                warning("In a supervised mode, the superblock corresponds to the response.")
        }
        return(FALSE)
    }else
        return(isTRUE(is_superblock))
}
check_tau <- function(tau, blocks, type = "rgcca",superblock=FALSE) {
    msg <- "tau should be comprise between 0 and 1 or should be set 'optimal' 
    for automatic setting"
    tau1 <- tau
    if(superblock){blocks[[length(blocks)+1]] <- Reduce(cbind,blocks);names(blocks)[length(blocks)]="superblock" }
    tryCatch({
        # Check value of each tau
        tau <- sapply(
            seq(length(tau)),
            function(x) {
                if (tau[x] != "optimal") {
                    y <- check_integer("tau", tau[x], float = TRUE, min = 0)
                    if (y > 1)
                        stop_rgcca(paste0(msg, " (currently equals to ", tau[x], ")."),
                                   exit_code = 129)
                    else
                        y
                }else
                    tau[x]
            })
        
        if (is(tau1, "matrix"))
            tau <- matrix(tau, NROW(tau1), NCOL(tau1))
        
        tau <- elongate_arg(tau, blocks)
        check_size_blocks(blocks, "tau", tau)
        tau <- check_spars(blocks, tau, type)
        
        return(tau)
        
        # If there is only one common tau
        # if (length(tau) == 1)
        #     tau <- rep(tau[[1]], length(blocks))
    }, warning = function(w)
        stop_rgcca(msg, exit_code = 131)
    )
}

