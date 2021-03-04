handle_NA <- function(blocks, method = "nipals") {
  if(method == "complete") blocks = intersection_list(blocks)
  if(method == "nipals")   blocks = blocks
  return(blocks)
}
