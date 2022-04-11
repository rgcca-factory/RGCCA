# Other checking functions
#-----------------------------
check_blockx <- function(x, y, blocks) {
  message <- paste0(
    x, " must be lower than the number of blocks, i.e. ",
    length(blocks), "."
  )
  exit_code <- 133
  x <- check_integer(x, y,
    max = length(blocks), exit_code = exit_code,
    max_message = message
  )
  return(x)
}

check_boolean <- function(x, y = x, type = "scalar") {
  if (is.null(y)) {
    y <- x
  }

  if (any(is.na(y))) {
    stop_rgcca(x, " must not be NA.")
  }

  if (!is(y, "logical")) {
    stop_rgcca(x, " must be TRUE or FALSE.")
  }

  if (type == "scalar" && length(y) != 1) {
    stop_rgcca(x, " must be of length 1.")
  }
}

check_colors <- function(colors) {
  if (!is.null(colors)) {
    colors <- as.vector(colors)
    lapply(colors, function(i) {
      if (!is.na(i) && !(i %in% colors()) && is.character2(i) &&
        regexpr("^#{1}[a-zA-Z0-9]{6,8}$", i) < 1) {
        stop_rgcca(
          "Unrecognized colors. Colors must be in colors() ",
          "or a rgb character."
        )
      }
    })
  }
}

check_compx <- function(x, y, ncomp, blockx) {
  res <- check_integer(x, y, min = 1)
  if (y > ncomp[blockx]) {
    stop_rgcca("not existing component. Trying to extract component ", y,
      " for block ", blockx, " , but only ", ncomp[blockx],
      " components are available for this block.",
      exit_code = 128
    )
  }
  return(res)
}

# Check the format of the connection matrix
check_connection <- function(C, blocks) {
  msg <- "connection matrix C must"

  if (!is.matrix(C)) stop_rgcca(msg, " be a matrix.", exit_code = 103)

  if (!isSymmetric.matrix(unname(C))) {
    stop_rgcca(paste(msg, "be symmetric."), exit_code = 103)
  }

  if (any(is.na(C))) {
    stop_rgcca(paste(msg, "not contain NA values."), exit_code = 106)
  }

  x <- C >= 0 & C <= 1
  if (sum(!x) != 0) {
    stop_rgcca(paste(msg, "contain numbers between 0 and 1."), exit_code = 106)
  }

  if (all(C == 0)) {
    stop_rgcca(paste(msg, "not contain only 0."), exit_code = 107)
  }

  invisible(check_size_blocks(blocks, "connection matrix", C))

  if (is.null(rownames(C)) || is.null(colnames(C))) {
    rownames(C) <- names(blocks) -> colnames(C)
  }

  if (!all(rownames(C) %in% names(blocks)) ||
    !all(colnames(C) %in% names(blocks))) {
    stop_rgcca(paste(
      msg,
      "have the rownames and the colnames that match with",
      "the names of the blocks."
    ),
    exit_code = 108
    )
  }

  return(C)
}

check_file <- function(f) {
  # Check the existence of a path f: A character giving the path of a file
  if (!file.exists(f)) {
    stop_rgcca(paste("file", f, "does not exist."), exit_code = 101)
  }
}

check_integer <- function(x, y = x, type = "scalar", float = FALSE, min = 1,
                          max = Inf, max_message = NULL, exit_code = NULL,
                          min_message = NULL) {
  if (is.null(y)) {
    y <- x
  }

  if (type %in% c("matrix", "data.frame")) {
    y_temp <- y
  }

  y <- tryCatch(
    as.double(as.matrix(y)),
    warning = function(w) {
      stop_rgcca(paste(x, "must be numeric."))
    }
  )

  if (any(is.na(y))) {
    stop_rgcca(paste(x, "must not be NA."))
  }

  if (type == "scalar" && length(y) != 1) {
    stop_rgcca(paste(x, "must be of length 1."))
  }

  if (!float) {
    if (any((y %% 1) != 0)) {
      stop_rgcca(paste(x, "must be an integer."))
    }
    y <- as.integer(y)
  }

  if (any(y < min)) {
    if (!is.null(min_message)) {
      stop_rgcca(min_message, exit_code = exit_code)
    } else {
      stop_rgcca(x, " must be higher than or equal to ", min, ".",
        exit_code = exit_code
      )
    }
  }

  if (any(y > max)) {
    if (!is.null(max_message)) {
      stop_rgcca(max_message, exit_code = exit_code)
    } else {
      stop_rgcca(x, " must be lower than or equal to ", max, ".",
        exit_code = exit_code
      )
    }
  }

  if (type %in% c("matrix", "data.frame")) {
    y <- matrix(
      y,
      NROW(y_temp),
      NCOL(y_temp),
      dimnames = dimnames(y_temp)
    )
  }

  if (type == "data.frame") {
    y <- as.data.frame(y)
  }

  return(y)
}

check_method <- function(method) {
  analysis <- c(
    "rgcca", "sgcca", "pca", "spca", "pls", "spls",
    "cca", "ifa", "ra", "gcca", "maxvar", "maxvar-b",
    "maxvar-a", "mcoa", "cpca-1", "cpca-2", "cpca-4",
    "hpca", "maxbet-b", "maxbet", "maxdiff-b", "maxdiff",
    "sabscor", "ssqcor", "ssqcov-1",
    "ssqcov-2", "ssqcov", "sumcor", "sumcov-1", "sumcov-2",
    "sumcov", "sabscov-1", "sabscov-2"
  )

  if (!tolower(method) %in% analysis) {
    stop_rgcca(
      "method '", method, "' is not among the available methods: ",
      paste(analysis, collapse = "', '"), "'.",
      exit_code = 112
    )
  }
  return(tolower(method))
}

check_nblocks <- function(blocks, method) {
  if (tolower(method) %in% c("pca", "spca")) {
    if (length(blocks) == 1) {
      return(blocks)
    }
    nb <- 1
    exit_code <- 110
  } else {
    if (length(blocks) == 2) {
      return(blocks)
    }
    nb <- 2
    exit_code <- 111
  }
  stop_rgcca(
    length(blocks),
    " blocks were provided but the number of blocks for ", method,
    " must be ", nb, ".",
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

check_ncomp <- function(ncomp, blocks, min = 1, superblock = FALSE,
                        response = NULL) {
  if (superblock) {
    if (length(unique(ncomp)) != 1) {
      stop_rgcca(
        "only one number of components must be specified (superblock)."
      )
    }
    max_ncomp <- ifelse(
      "superblock" %in% names(blocks),
      NCOL(blocks[[length(blocks)]]),
      sum(vapply(blocks, NCOL, FUN.VALUE = integer(1)))
    )
    msg <- paste0(
      "the number of components must be lower than the number of ",
      "variables in the superblock, i.e. ", max_ncomp,
      "."
    )

    y <- check_integer("ncomp", ncomp[1],
      min = min, max_message = msg,
      max = max_ncomp,
      exit_code = 126
    )
    return(rep(y, length(ncomp)))
  }

  ncomp <- elongate_arg(ncomp, blocks)
  check_size_blocks(blocks, "ncomp", ncomp)
  ncomp <- vapply(
    seq(length(ncomp)),
    function(x) {
      if (!is.null(response) && x == response) {
        y <- check_integer("ncomp", ncomp[x], min = min, exit_code = 126)
      } else {
        msg <- paste0(
          "ncomp[", x, "] must be lower than the number of variables ",
          "for block ", x, ", i.e. ", NCOL(blocks[[x]]), "."
        )
        y <- check_integer("ncomp", ncomp[x],
          min = min, max_message = msg,
          max = NCOL(blocks[[x]]), exit_code = 126
        )
      }
      return(y)
    },
    FUN.VALUE = integer(1)
  )
  return(ncomp)
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
      "contains qualitative data. Transform them in a disjunctive table."
    )

    if (!header) {
      msg <- paste0(
        msg, "Possible mistake: header parameter is disabled, ",
        "check if the file doesn't have one."
      )
    }

    stop_rgcca(paste(msg, "\n"), exit_code = 100)
  }
}

check_response <- function(response = NULL) {
  if (is.null(response)) {
    return(NA)
  }
  response <- as.data.frame(response)
  qualitative <- is.character2(response)
  if (qualitative || (NCOL(response) == 1)) {
    return(response)
  }

  # If response is numeric, contains only one nonzero value per row that equals
  # 1, we convert response to a column factor.
  response <- to_numeric(response)
  col_sum <- unique(apply(response, 1, sum))
  values <- unique(unlist(response))
  values <- union(col_sum, values)
  if (identical(sort(union(values, c(0, 1))), c(0, 1))) {
    response <- factor(
      apply(response, 1, which.max),
      levels = colnames(response)
    )
    return(
      data.frame(
        as.character(response),
        row.names = rownames(response)
      )
    )
  }
  return(response)
}

# Test on the sign of the correlation
check_sign_comp <- function(rgcca_res, w) {
  y <- lapply(
    seq_along(rgcca_res$a),
    function(i) pm(rgcca_res$call$blocks[[i]], w[[i]])
  )

  w <- lapply(setNames(seq_along(w), names(w)), function(i) {
    if (NROW(w[[i]]) < NROW(y[[i]])) {
      res <- as.matrix(cor(rgcca_res$Y[[i]], y[[i]]))
    } else {
      res <- as.matrix(cor(rgcca_res$a[[i]], w[[i]]))
    }
    vec_sign <- vapply(diag(res), function(x) {
      return(ifelse(!is.na(x) && (x < 0), -1, 1))
    }, double(1))
    return(pm(w[[i]], diag(vec_sign, nrow = nrow(res))))
  })

  return(w)
}

check_size_blocks <- function(blocks, x, y = x, n_row = NULL) {
  if (identical(x, y)) {
    x <- ""
  }
  if (any(class(y) %in% c("matrix", "data.frame"))) {
    dim_y <- NCOL(y)
    dim_type <- "number of columns"
    if (!is.null(n_row) && (NROW(y) != n_row) && (NROW(y) != 1)) {
      stop_rgcca(x, " must have ", n_row, " rows.")
    }
  } else {
    dim_y <- length(y)
    dim_type <- "size"
  }

  if (dim_y != length(blocks)) {
    stop_rgcca(
      x,
      " must have the same ",
      dim_type,
      " (actually ",
      dim_y,
      ") as the number of blocks (",
      length(blocks),
      ").",
      exit_code = 130
    )
  } else {
    return(TRUE)
  }
}


# Print warning if file size over
check_size_file <- function(filename) {
  size <- file.size(filename)
  if (size > 5e+06) {
    # warning(paste0('The size of ', filename, ' is over 5 Mo (',
    #  round(size / 1E6, 1), ' Mo). File loading could take some times...'),
    message("File loading in progress ...")
  }
}

check_penalty <- function(penalty, blocks, method = "rgcca", superblock = FALSE,
                          ncomp = NULL) {
  if (superblock) {
    blocks[[length(blocks) + 1]] <- Reduce(cbind, blocks)
    names(blocks)[length(blocks)] <- "superblock"
  }
  penalty <- elongate_arg(penalty, blocks)
  name <- ifelse(method == "rgcca", "tau", "sparsity")
  check_size_blocks(blocks, name, penalty, n_row = ncomp)

  is_matrix <- is.matrix(penalty)
  DIM <- dim(penalty)

  # Check value of each penalty
  if (method == "rgcca") {
    penalty <- unlist(lapply(penalty, check_tau))
  }
  if (method == "sgcca") {
    divider <- ifelse(is_matrix, DIM[1], 1)
    penalty <- vapply(
      seq(length(penalty)),
      function(x) {
        n <- 1 + (x - 1) / divider
        check_spars(penalty[x], blocks[[n]], n)
      }, FUN.VALUE = double(1L)
    )
  }

  if (is_matrix) penalty <- matrix(penalty, DIM[1], DIM[2])

  return(penalty)
}

check_spars <- function(sparsity, block, n) {
  min_sparsity <- 1 / sqrt(NCOL(block))
  min_message <- paste0(
    "too high sparsity. Sparsity parameter equals ", sparsity,
    ". For SGCCA, it must be greater than ",
    "1/sqrt(number_column) (i.e., ", round(min_sparsity, 4),
    " for block ", n, ")."
  )
  sparsity <- check_integer("sparsity", sparsity,
    float = TRUE,
    min = min_sparsity, max = 1, min_message = min_message
  )
  invisible(sparsity)
}

check_tau <- function(tau) {
  if (is.na(tau) || tau != "optimal") {
    tau <- check_integer("tau", tau, float = TRUE, min = 0, max = 1)
  }
  invisible(tau)
}

check_scheme <- function(scheme) {
  if (
    (mode(scheme) != "function") &&
      (scheme != "horst") &&
      (scheme != "factorial") &&
      (scheme != "centroid")
  ) {
    stop_rgcca(paste0(
      "scheme must be one of the following schemes: horst, ",
      "centroid, factorial or a function."
    ))
  }
}
