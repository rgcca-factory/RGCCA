available_NA_methods <- c("complete", "nipals")

handle_NA <- function(blocks, NA_method = "nipals") {
  na.rm <- FALSE
  if (!NA_method %in% available_NA_methods) {
    stop_rgcca(paste0(
      "NA_method ", NA_method, " is not implemented to handle missing values.",
      "Please select one among (",
      paste(available_NA_methods, collapse = ", "), ")."
    ))
  }
  if (NA_method == "complete") blocks <- intersection_list(blocks)
  if (NA_method == "nipals") {
    blocks <- blocks
    na.rm <- Reduce("||", lapply(blocks, function(x) any(is.na(x))))
  }
  if(NA_method == "imputeSB")
  {
    blocks=impute_sb(blocks,
                     tau = rep(1, length(blocks) + 1),
                     ni = 5,
                     tol = 1e-08,
                     ncomp = NULL,
                     naxis = 1,
                     scale = TRUE,
                     scale_block = TRUE,
                     scheme = "centroid",
                     bias = TRUE,
                     verbose = FALSE,
                     threshold = 1e-6,
                     reg = "y",
                     quiet = FALSE)
  }
  if(NA_method == "impute")
  {
    blocks=impute(blocks,
                    connection,
                    tau = rep(1, length(blocks) + 1),
                    ni = 5,
                    tol = 1e-08,
                    ncomp = NULL,
                    naxis = 1,
                    scale = TRUE,
                    scale_block = TRUE,
                    scheme = "centroid",
                    bias = TRUE,
                    verbose = FALSE,
                    threshold = 1e-6,
                    reg = "y",
                    quiet = FALSE)
  }
  return(list(blocks = blocks, na.rm = na.rm))
}
