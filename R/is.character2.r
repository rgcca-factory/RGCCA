# Test if at least a character is present in a collection of elements.
#
# Tests if a collection of elements is composed only by qualitative variables.
# To do so, we try to cast the elements to numeric values. If it does not work,
# NA values are produced so we return TRUE in this case.
#
# @param x A collection to be tested.
# @return A bolean for the presence (FALSE) or the absence (TRUE) of at least
# one quantitative variable.
# @param warn_separator A boolean to print warning for bad separator use.
is.character2 <- function(x, warn_separator = FALSE) {
  # is. character() consider a string with '1.2' as a character, not this
  # function.
  local_env <- new.env()
  assign("res", FALSE, envir = local_env)
  x <- unlist(x)
  x <- x[x != "NA"]
  tryCatch(
    as.numeric(x),
    warning = function(w) assign("res", TRUE, envir = local_env)
  )

  res <- get("res", envir = local_env)

  if (warn_separator && res) {
    stop_rgcca("Wrong separator. Please, select the comma separator parameter.")
  }

  return(res)
}
