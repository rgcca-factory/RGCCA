#' Small utility function to format a matrix of parameters
#' @noRd
format_combinations <- function(par_value) {
  
  combinations=list()
 
  for (i in 1:NROW(par_value)){
    combinations[i] = paste0(round(unlist(par_value[i, ]), 2),collapse='/')
  }
  combinations=unlist(combinations)

  # If parameters are too long, there are replaced with "Set x"
  # The same is done if rounding to 2 digits leads to the same values
  to_set <- (nchar(combinations[1]) > 15) |
    (length(unique(combinations)) < NROW(par_value))
  if (to_set) {
    combinations <- paste("Set ", sep = "", seq_len(NROW(par_value)))
  }
  return(combinations)
}
