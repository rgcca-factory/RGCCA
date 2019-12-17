scaling = function(
    blocks,
    scale = TRUE,
    bias = TRUE) {

    if (scale) {
        lapply(blocks, function(x)
        scale2(x, bias = bias) / sqrt(NCOL(x)))
    }else{
        lapply(blocks, function(x)
        scale2(x, scale = FALSE))
    }
}
