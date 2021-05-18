# Keeps only subject without missing values
# @inheritParams rgccad
# @return The intersected list from matrices
# @title intersection_list
# @examples
# set.seed(42);X1=matrix(rnorm(35),7,5);
# set.seed(22);X2=matrix(rnorm(28),7,4);
# set.seed(2);X3=matrix(rnorm(49),7,7);
## usual test
#X1[1,]=NA
#X2[7,1]=NA
#X2[5,1]=NA
#A=list(X1,X2)
#intersection_list(A=A)
## too many subjects with missing values
#X3[3,1:2]=NA
#intersection_list(A=list(X1,X2,X3))
# @export intersection_list
intersection_list <- function(A) {
  # Find rows without missing values in each block
  valid_rows        = lapply(A, complete.cases)
  # Take the intersection
  common_valid_rows = apply(
    matrix(unlist(valid_rows), length(valid_rows[[1]]), length(valid_rows)),
    1, prod)
  if (sum(common_valid_rows) <= 3) {
    stop_rgcca(paste0("Less than 3 subjects have no missing values, choose",
                      " another missing value handling method or work on ",
                      "your dataset."))
  }
  # Extract the rows from the different blocks
  lapply(A, function(x) subset_rows(x, as.logical(common_valid_rows)))
}
