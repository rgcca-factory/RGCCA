# File name from a path
#
# Get the file name from a path
#
# @param fi A character giving the path of a file
# @return A character for the name of the file

get_filename <- function(fi) {
  if (!is.null(fi)) {
    fo <- unlist(strsplit(fi, "/"))
    fo <- fo[length(fo)]
    unlist(strsplit(fo, "[.]"))[1]
  }
}
