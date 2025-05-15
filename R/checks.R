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

check_colors <- function(colors, type = "variables") {
  if (is.null(colors)) {
    switch(type,
      "variables" = colors <- c(
        "#882255", "#1FB0B5", "#C79E46", "#117733",
        "#332288", "#BB728E", "#1F9F76", "#A95700"
      ),
      "samples" = colors <- c(
        "#1F6EE0", "#CC0000", "#029B38", "#FFC709",
        "#34C926", "#88AEC1", "#990099", "#FB176C"
      ),
      "AVE" = colors <- c(
        "#828076", "#959685", "#A4AC96", "#AEB998",
        "#B7C69A", "#CADF9E", "#DCF1AE", "#E7F8C4"
      )
    )
  } else {
    colors <- as.vector(colors)
    lapply(colors, function(i) {
      if (!is.na(i) && !(i %in% colors()) && is.character(i) &&
        regexpr("^#{1}[a-zA-Z0-9]{6,8}$", i) < 1) {
        stop_rgcca(
          "Unrecognized colors. Colors must be in colors() ",
          "or a rgb character."
        )
      }
    })
  }
  return(colors)
}

check_compx <- function(x, y, ncomp, blockx) {
  res <- check_integer(x, y, min = 1)
  if (y > ncomp[blockx]) {
    stop_rgcca("not existing component. Trying to extract component ", y,
      " for block ", blockx, " , but only ", ncomp[blockx],
      " component(s) are available for this block.",
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
    rownames(C) <- colnames(C) <- names(blocks)
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

check_integer <- function(x, y = x, type = "scalar", float = FALSE, min = 1,
                          max = Inf, max_message = NULL, exit_code = NULL,
                          min_message = NULL) {
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
  if (!method %in% available_methods()) {
    stop_rgcca(
      "method '", method, "' is not among the available methods: ",
      paste(available_methods(), collapse = "', '"), "'.",
      exit_code = 112
    )
  }
  return(method)
}

check_nblocks <- function(blocks, method) {
  if (tolower(method) %in% one_block_methods()) {
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
    seq_along(ncomp),
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

# Test on the sign of the correlation
check_sign_comp <- function(rgcca_res, w) {
  y <- lapply(
    seq_along(rgcca_res$a),
    function(i) pm(rgcca_res$blocks[[i]], w[[i]])
  )

  w <- lapply(setNames(seq_along(w), names(w)), function(i) {
    if (NROW(w[[i]]) < NROW(y[[i]])) {
      res <- as.matrix(cor2(rgcca_res$Y[[i]], y[[i]]))
    } else {
      res <- as.matrix(cor2(rgcca_res$a[[i]], w[[i]]))
    }
    vec_sign <- vapply(diag(res), function(x) {
      return(ifelse(!is.na(x) && (x < 0), -1, 1))
    }, double(1))
    return(pm(w[[i]], diag(vec_sign, nrow = nrow(res))))
  })

  return(w)
}

check_size_blocks <- function(blocks, x, y = x, n_row = NULL,
                              superblock = FALSE) {
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
    superblock_msg <- ifelse(superblock, paste0(
      " or the number of blocks + 1 (", length(blocks) + 1, ")"
    ), "")
    message <- paste0(
      x,
      " must have the same ",
      dim_type,
      " (actually ",
      dim_y,
      ") as the number of blocks (",
      length(blocks),
      ")", superblock_msg, "."
    )
    stop_rgcca(message, exit_code = 130
    )
  } else {
    return(TRUE)
  }
}

check_penalty <- function(penalty, blocks, method = "rgcca", superblock = FALSE,
                          ncomp = NULL) {
  penalty <- elongate_arg(penalty, blocks)
  is_matrix <- is.matrix(penalty)
  DIM <- dim(penalty)
  size <- ifelse(is_matrix, NCOL(penalty), NROW(penalty))
  if (superblock && (size == (length(blocks) + 1))) {
    blocks[[length(blocks) + 1]] <- Reduce(cbind, blocks)
    names(blocks)[length(blocks)] <- "superblock"
  }
  name <- ifelse(method == "rgcca", "tau", "sparsity")
  check_size_blocks(blocks, name, penalty,
                    n_row = ncomp, superblock = superblock)

  # Check value of each penalty
  if (method == "rgcca") {
    penalty <- unlist(lapply(penalty, check_tau))
  }
  if (method == "sgcca") {
    divider <- ifelse(is_matrix, DIM[1], 1)
    penalty <- vapply(
      seq_along(penalty),
      function(x) {
        n <- 1 + (x - 1) / divider
        check_spars(penalty[x], blocks[[n]], n)
      },
      FUN.VALUE = double(1L)
    )
  }

  if (is_matrix) penalty <- matrix(penalty, DIM[1], DIM[2])

  return(penalty)
}

check_spars <- function(sparsity, block, n) {
  if (mode(block) == "character") {
    return(sparsity)
  }
  min_sparsity <- 1 / sqrt(NCOL(block))
  min_message <- paste0(
    "too low sparsity. Sparsity parameter equals ", sparsity,
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
  if (mode(scheme) != "function") {
    scheme <- tolower(scheme)
    if (!scheme %in% c("horst", "factorial", "centroid")) {
      stop_rgcca(paste0(
        "scheme must be one of the following schemes: 'horst', ",
        "'centroid', 'factorial' or a function."
      ))
    }
  }
  return(scheme)
}

check_prediction_model <- function(prediction_model, response_block,
                                   missing_model = FALSE) {
  classification <-
    is.factor(response_block) || is.character(response_block)

  if (missing_model && classification) {
    prediction_model <- "lda"
  }

  if (is.list(prediction_model)) {
    model_info <- prediction_model
  } else {
    model_info <- caret::getModelInfo(prediction_model, regex = FALSE)[[1]]
    if (is.null(model_info)) {
      stop_rgcca(
        "unknown model. Model ", prediction_model, " is not handled, please ",
        "see caret::modelLookup() for a list of the available models."
      )
    }
  }

  is_inadequate <- !("Classification" %in% model_info$type) && classification
  if (is_inadequate) {
    stop_rgcca(
      "inadequate model. Response block contains categorical data ",
      "but model ", prediction_model, " is not made for ",
      "classification. Please choose another model."
    )
  }

  is_inadequate <- !("Regression" %in% model_info$type) && !classification
  if (is_inadequate) {
    stop_rgcca(
      "inadequate model. Response block contains continuous data ",
      "but model ", prediction_model, " is not made for ",
      "regression Please choose another model."
    )
  }

  return(list(prediction_model = model_info, classification = classification,
              model_name = prediction_model))
}

check_char <- function(arg, name_arg, values) {
  res <- grep(arg, values, fixed = TRUE, value = TRUE)
  if (length(res) == 0) {
    stop_rgcca(
      "'", name_arg, "' should be one of \"",
      paste(values, collapse = "\", \""), "\""
    )
  } else {
    return(res[1])
  }
}

is.sparseMatrix <- function(x) is(x, 'sparseMatrix')
check_matrix <- function(x) is.matrix(x) || is.sparseMatrix(x)

check_laplacians <- function(laplacians, blocks) {
  if (is.null(laplacians) || all(sapply(laplacians, is.null))) {
    return(NULL)
  }
  if (!is.list(laplacians) || (length(laplacians)!=length(blocks))) {
      stop_rgcca("graph_laplacians should be a list",
        "with as many elements as blocks")
  }
  lap_not_null = laplacians[!sapply(laplacians, is.null)]
  if (!all(sapply(lap_not_null, check_matrix))) {
    stop_rgcca("all non NULL elements of graph_laplacians",
      "should be either regular or sparse matrices")
  }
  return(laplacians)
}
